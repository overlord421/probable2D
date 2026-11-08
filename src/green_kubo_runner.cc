#include <green_kubo_runner.hh>
#include <pauli_time_histogram.hh>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <vector>

#if __has_include(<omp.h>)
#include <omp.h>
#endif

namespace probable {

namespace {

void wrap_periodic(const Material &material, Pos2D &r) {
  while (r.x < 0) r.x += material.Lx;
  while (r.x >= material.Lx) r.x -= material.Lx;
  while (r.y < 0) r.y += material.Ly;
  while (r.y >= material.Ly) r.y -= material.Ly;
}

double axis_field(const Material &material, const Vec2 &field) {
  return material.thermal_axis == 'y' ? field.y : field.x;
}

double stddev(const std::vector<double> &values) {
  if (values.size() < 2) {
    return 0;
  }
  double mean = 0;
  double m2 = 0;
  size_t n = 0;
  for (double value : values) {
    ++n;
    double delta = value - mean;
    mean += delta / n;
    m2 += delta * (value - mean);
  }
  return std::sqrt(m2 / (values.size() - 1));
}

struct Correlations {
  std::vector<double> qq;
  std::vector<double> qj;
  std::vector<double> jj;
};

std::unique_ptr<PauliTimeHistogram> prepare_pauli_histogram(
    const Material &material, const std::vector<Scattering *> &mechanisms,
    const GreenKuboConfig &config) {
  if (!config.pauli_time_histogram) return nullptr;
  if (!material.use_fermi_dirac || config.electric_field.x != 0 ||
      config.electric_field.y != 0 || config.magnetic_field_z != 0)
    throw std::invalid_argument("Time-histogram Pauli mode requires FD and zero fields");
  if (!config.pauli_preparation_passes || !config.pauli_histories_per_pass ||
      !(config.time_step > 0) || !(config.temperature > 0) ||
      !(config.pauli_preparation_time > 0) || !std::isfinite(config.pauli_preparation_time) ||
      std::fabs(config.carrier_density_2d / material.carrier_density_2d - 1) > 1e-12)
    throw std::invalid_argument("Invalid Pauli preparation configuration");
  double excess = std::max(0.0, material.chemical_potential(config.temperature) - material.Delta)
      + 12 * consts::kB * config.temperature;
  double px_limit = std::sqrt(material.mx * excess * (2 * material.Delta + excess) / material.Delta);
  double py_limit = std::sqrt(2 * material.my * excess);
  auto histogram = std::make_unique<PauliTimeHistogram>(config.pauli_bins, px_limit, py_limit,
      config.carrier_density_2d, material.spin_degeneracy);
  std::vector<double> reference(histogram->cells());
  for (size_t i = 0; i < histogram->cells(); ++i) {
    Vec2 center = histogram->center(i);
    double sum = 0;
    for (size_t x = 0; x < 16; ++x) for (size_t y = 0; y < 16; ++y) {
      Vec2 p{center.x + ((x + .5) / 16 - .5) * histogram->dp_x(),
             center.y + ((y + .5) / 16 - .5) * histogram->dp_y()};
      sum += material.occupation(material.energy(p), config.temperature);
    }
    reference[i] = sum / 256;
    histogram->initialize(i, reference[i]);
  }
  if (!config.pauli_occupation_file.empty()) {
    std::ifstream input(config.pauli_occupation_file);
    if (!input) throw std::invalid_argument("Cannot open Pauli occupation map");
    std::vector<double> raw;
    double inferred_total = 0, largest_raw = 0;
    std::string line;
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '#') continue;
      size_t i = raw.size();
      double px, py, residence, f_raw, f_used, fd;
      std::istringstream row(line);
      if (!(row >> px >> py >> residence >> f_raw >> f_used >> fd) || i >= histogram->cells())
        throw std::invalid_argument("Malformed Pauli occupation map");
      Vec2 center = histogram->center(i);
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(f_raw) || !std::isfinite(fd) ||
          std::fabs(px - center.x) > px_limit * 1e-10 ||
          std::fabs(py - center.y) > py_limit * 1e-10 || std::fabs(fd - reference[i]) > 1e-8)
        throw std::invalid_argument("Pauli map geometry or material parameters differ from run");
      histogram->initialize(i, f_used);
      if (residence < 0 || !std::isfinite(residence))
        throw std::invalid_argument("Invalid residence time in Pauli map");
      if (residence > 0) histogram->record(center, residence * units::s);
      if (f_raw > largest_raw && residence > 0) {
        largest_raw = f_raw;
        inferred_total = config.carrier_density_2d * residence * units::s
            / histogram->states_per_area() / f_raw;
      }
      raw.push_back(f_raw);
    }
    if (raw.size() != histogram->cells() || !(histogram->total_time() > 0))
      throw std::invalid_argument("Incomplete Pauli occupation map");
    if (!(inferred_total > 0) || inferred_total < histogram->total_time() * (1 - 1e-8))
      throw std::invalid_argument("Inconsistent total residence in Pauli map");
    if (inferred_total > histogram->total_time())
      histogram->record({2 * px_limit, 0}, inferred_total - histogram->total_time());
    if (histogram->outside_fraction() > 1e-3)
      throw std::invalid_argument("Pauli map has unresolved out-of-grid visits");
    for (size_t i = 0; i < raw.size(); ++i) {
      if (std::fabs(raw[i] - histogram->raw_occupation(i)) > 1e-8)
        throw std::invalid_argument("Pauli map density, state capacity or overflow is inconsistent");
      Vec2 p = histogram->center(i);
      if (std::fabs(histogram->occupation(p) - histogram->occupation({-p.x, -p.y})) > 1e-10 ||
          std::fabs(histogram->occupation(p) - histogram->occupation({-p.x, p.y})) > 1e-10)
        throw std::invalid_argument("Pauli map breaks zero-field reflection symmetry");
    }
    if (std::fabs(histogram->clipped_density_fraction() - 1) > .01)
      throw std::invalid_argument("Pauli map density normalization is unresolved");
    std::filesystem::create_directories(config.output_dir);
    std::string destination = config.output_dir + "/pauli_occupation_"
        + (material.thermal_axis == 'y' ? "y" : "x") + ".tsv";
    if (std::filesystem::absolute(config.pauli_occupation_file).lexically_normal() !=
        std::filesystem::absolute(destination).lexically_normal())
      std::filesystem::copy_file(config.pauli_occupation_file, destination,
                                std::filesystem::copy_options::overwrite_existing);
    return histogram;
  }
  size_t steps = static_cast<size_t>(config.pauli_preparation_time / config.time_step /
      config.pauli_preparation_passes / config.pauli_histories_per_pass);
  size_t stride = std::max<size_t>(1, config.sample_stride);
  if (steps < 10 * stride)
    throw std::invalid_argument("Pauli preparation histories are too short for time sampling");
  std::filesystem::create_directories(config.output_dir);
  std::string suffix = material.thermal_axis == 'y' ? "y" : "x";
  std::ofstream passes(config.output_dir + "/pauli_preparation_" + suffix + ".tsv");
  passes << std::setprecision(17)
      << "# pass cumulative_residence_s max_map_change outside_fraction used_density_over_n raw_saturated_cells\n";
  // Successive histories are run one at a time. Independent starting energies
  // cover the conserved energy modulo the single optical phonon quantum.
  // FD initializes particles and the first pass; later blocking uses counts.
  double last_change = 0;
  for (size_t pass = 0; pass < config.pauli_preparation_passes; ++pass) {
    for (size_t history = 0; history < config.pauli_histories_per_pass; ++history) {
      Vec2 p = material.create_particle_at_temperature(config.temperature);
      std::vector<double> free_flight(mechanisms.size());
      for (double &l : free_flight) l = -std::log(uniform());
      // During zero-field preparation p and the sampling coordinate stay
      // unchanged until an accepted collision. Cache those identical rates;
      // keep the original time-step hazards and random-event ordering.
      std::vector<double> preparation_rates(mechanisms.size());
      for (size_t k = 0; k < mechanisms.size(); ++k)
        preparation_rates[k] = mechanisms[k]->rate(p, 0);
      size_t burn = steps / 5;
      for (size_t step = 0; step < steps; ++step) {
        if (step >= burn && (step - burn) % stride == 0)
          histogram->record_reflected(p, std::min(stride, steps - step) * config.time_step);
        for (size_t k = 0; k < mechanisms.size(); ++k) {
          free_flight[k] -= preparation_rates[k] * config.time_step;
          if (free_flight[k] < 0) {
            Vec2 final = mechanisms[k]->scatter(p);
            double f = pass == 0 ? material.occupation(material.energy(final), config.temperature)
                                 : histogram->occupation(final);
            if (uniform() >= f) {
              p = final;
              for (size_t j = 0; j < mechanisms.size(); ++j)
                preparation_rates[j] = mechanisms[j]->rate(p, 0);
            }
            free_flight[k] = -std::log(uniform());
            break;
          }
        }
      }
    }
    // Under-relaxation prevents noisy occupancy estimates from turning a
    // nearly full cell into an artificial absorbing region. Scaling by the
    // previous vacancy stabilizes the strongly degenerate fixed point.
    last_change = histogram->update(.5, true);
    passes << pass << ' ' << histogram->total_time() / units::s << ' ' << last_change
           << ' ' << histogram->outside_fraction() << ' ' << histogram->clipped_density_fraction()
           << ' ' << histogram->saturated_cells() << '\n';
  }
  std::ofstream map(config.output_dir + "/pauli_occupation_" + suffix + ".tsv");
  map << std::setprecision(17)
      << "# px_internal py_internal residence_s f_raw f_used f_FD_cell_average\n";
  for (size_t i = 0; i < histogram->cells(); ++i) {
    Vec2 p = histogram->center(i);
    map << p.x << ' ' << p.y << ' ' << histogram->residence(i) / units::s << ' '
        << histogram->raw_occupation(i) << ' ' << histogram->frozen_occupation(i) << ' '
        << reference[i] << '\n';
  }
  map.close();
  passes.close();
  if (histogram->outside_fraction() > 1e-3 ||
      std::fabs(histogram->clipped_density_fraction() - 1) > .01 || last_change > .02)
    throw std::runtime_error("Pauli preparation is unresolved: inspect pauli_preparation_*.tsv; "
                            "increase preparation time/passes or coarsen the momentum grid");
  return histogram;
}

Correlations compute_correlations(const std::vector<double> &heat_current,
                                  const std::vector<double> &particle_current,
                                  size_t ensemble_size,
                                  size_t max_lag) {
  size_t samples = heat_current.size();
  size_t lags = std::min(max_lag, samples > 0 ? samples - 1 : 0);
  Correlations corr;
  corr.qq.assign(lags + 1, 0);
  corr.qj.assign(lags + 1, 0);
  corr.jj.assign(lags + 1, 0);
  if (samples == 0 || ensemble_size == 0) {
    return corr;
  }

  double inv_ensemble = 1.0 / ensemble_size;
  for (size_t lag = 0; lag <= lags; ++lag) {
    double qq = 0;
    double qj = 0;
    double jj = 0;
    size_t count = samples - lag;
    for (size_t i = 0; i < count; ++i) {
      double q0 = heat_current[i] * inv_ensemble;
      double q1 = heat_current[i + lag] * inv_ensemble;
      double j0 = particle_current[i] * inv_ensemble;
      double j1 = particle_current[i + lag] * inv_ensemble;
      qq += q1 * q0;
      qj += q1 * j0;
      jj += j1 * j0;
    }
    double scale = static_cast<double>(ensemble_size) / count;
    corr.qq[lag] = qq * scale;
    corr.qj[lag] = qj * scale;
    corr.jj[lag] = jj * scale;
  }
  return corr;
}

struct RunningKappa {
  std::vector<double> iqq;
  std::vector<double> iqj;
  std::vector<double> ijj;
  std::vector<double> kappa;
  std::vector<double> kappa_j0;
};

RunningKappa integrate_correlations(const Correlations &corr,
                                    double sample_dt,
                                    double temperature,
                                    double carrier_density_2d) {
  size_t lags = corr.qq.size();
  RunningKappa running;
  running.iqq.assign(lags, 0);
  running.iqj.assign(lags, 0);
  running.ijj.assign(lags, 0);
  running.kappa.assign(lags, 0);
  running.kappa_j0.assign(lags, 0);
  if (lags == 0) {
    return running;
  }

  double prefactor = carrier_density_2d / (consts::kB * temperature * temperature);
  double kappa_unit = units::J / units::s / units::K;
  for (size_t lag = 0; lag < lags; ++lag) {
    double weight = lag == 0 ? 0.5 : 1.0;
    running.iqq[lag] = (lag == 0 ? 0 : running.iqq[lag - 1]) + weight * corr.qq[lag] * sample_dt;
    running.iqj[lag] = (lag == 0 ? 0 : running.iqj[lag - 1]) + weight * corr.qj[lag] * sample_dt;
    running.ijj[lag] = (lag == 0 ? 0 : running.ijj[lag - 1]) + weight * corr.jj[lag] * sample_dt;

    double corrected = running.iqq[lag];
    if (std::fabs(running.ijj[lag]) > 1e-300) {
      corrected -= running.iqj[lag] * running.iqj[lag] / running.ijj[lag];
    }
    running.kappa[lag] = prefactor * running.iqq[lag] / kappa_unit;
    running.kappa_j0[lag] = prefactor * corrected / kappa_unit;
  }
  return running;
}

} // namespace

GreenKuboRunResult run_green_kubo_simulation(const Material &material,
                                             const std::vector<Scattering *> &mechanisms,
                                             const GreenKuboConfig &config) {
  auto pauli_histogram = prepare_pauli_histogram(material, mechanisms, config);
  size_t steps = config.all_time / config.time_step + 1;
  size_t tally_start_step = steps / 5;
  size_t effective_steps = steps > tally_start_step ? steps - tally_start_step : 0;
  size_t sample_stride = std::max<size_t>(1, config.sample_stride);
  size_t samples = (effective_steps + sample_stride - 1) / sample_stride;
  double sample_dt = sample_stride * config.time_step;

  std::vector<double> heat_current(samples, 0);
  std::vector<double> particle_current(samples, 0);
  std::vector<double> scattering_counts(mechanisms.size(), 0);

#pragma omp parallel
  {
    std::vector<double> local_heat(samples, 0);
    std::vector<double> local_particle(samples, 0);
    std::vector<double> local_scattering_counts(mechanisms.size(), 0);

#pragma omp for
    for (size_t i = 0; i < config.ensemble_size; ++i) {
      Pos2D r = material.create_initial_position();
      Vec2 p = material.create_particle_at_temperature(config.temperature);
      std::vector<double> free_flight(mechanisms.size(), 0);
      for (double &l : free_flight) {
        l = -std::log(uniform());
      }

      for (size_t step = 0; step < steps; ++step) {
        Vec2 v = material.velocity(p);
        double e = material.energy(p);

        if (step >= tally_start_step && ((step - tally_start_step) % sample_stride == 0)) {
          size_t sample = (step - tally_start_step) / sample_stride;
          double v_axis = material.velocity_component(v);
          local_heat[sample] += (e - material.Delta) * v_axis;
          local_particle[sample] += v_axis;
        }

        r.x += v.x * config.time_step;
        r.y += v.y * config.time_step;
        wrap_periodic(material, r);

        p += -consts::e * (config.electric_field + v.cross_with_B(config.magnetic_field_z)) * config.time_step;

        for (size_t k = 0; k < mechanisms.size(); ++k) {
          free_flight[k] -= mechanisms[k]->rate(p, material.axis_coordinate(r)) * config.time_step;
          if (free_flight[k] < 0) {
            Vec2 scattered = mechanisms[k]->scatter(p);
            double final_occupation =
                !material.use_fermi_dirac ? 0 : pauli_histogram ? pauli_histogram->occupation(scattered)
                : material.occupation(material.energy(scattered), config.temperature);
            bool accepted = !material.use_fermi_dirac || uniform() >= final_occupation;
            if (accepted) {
              p = scattered;
              local_scattering_counts[k] += 1;
            }
            free_flight[k] = -std::log(uniform());
            break;
          }
        }
      }
    }

#pragma omp critical
    {
      for (size_t i = 0; i < samples; ++i) {
        heat_current[i] += local_heat[i];
        particle_current[i] += local_particle[i];
      }
      for (size_t i = 0; i < mechanisms.size(); ++i) {
        scattering_counts[i] += local_scattering_counts[i];
      }
    }
  }

  Correlations corr =
      compute_correlations(heat_current, particle_current, config.ensemble_size, config.max_correlation_lag);
  RunningKappa running =
      integrate_correlations(corr, sample_dt, config.temperature, config.carrier_density_2d);

  GreenKuboRunResult result;
  if (pauli_histogram) {
    result.pauli_histogram_outside_fraction = pauli_histogram->outside_fraction();
    result.pauli_histogram_clipped_density_fraction = pauli_histogram->clipped_density_fraction();
    result.pauli_histogram_saturated_cells = pauli_histogram->saturated_cells();
  }
  result.samples = samples;
  result.correlation_lags = corr.qq.size();
  result.scattering_counts = scattering_counts;
  result.scattering_rates.assign(mechanisms.size(), 0);
  for (size_t i = 0; i < mechanisms.size(); ++i) {
    result.scattering_rates[i] = scattering_counts[i] / config.ensemble_size / config.all_time * units::s;
  }

  double heat_flux_unit_2d = units::J / units::s / units::m;
  double avg_heat = 0;
  double avg_particle = 0;
  for (size_t i = 0; i < samples; ++i) {
    avg_heat += heat_current[i] / config.ensemble_size;
    avg_particle += particle_current[i] / config.ensemble_size;
  }
  if (samples > 0) {
    avg_heat /= samples;
    avg_particle /= samples;
  }
  result.heat_flux_2d = avg_heat * config.carrier_density_2d / heat_flux_unit_2d;
  result.particle_flux_2d = avg_particle * config.carrier_density_2d;
  if (!running.kappa.empty()) {
    result.kappa_2d = running.kappa.back();
    result.kappa_j0_2d = running.kappa_j0.back();
  }

  size_t blocks = std::min(config.blocks, samples);
  std::vector<double> block_kappa;
  std::vector<double> block_kappa_j0;
  if (blocks > 1) {
    size_t block_size = samples / blocks;
    for (size_t block = 0; block < blocks; ++block) {
      size_t begin = block * block_size;
      size_t end = block == blocks - 1 ? samples : begin + block_size;
      if (end <= begin + 1) {
        continue;
      }
      std::vector<double> block_heat(heat_current.begin() + begin, heat_current.begin() + end);
      std::vector<double> block_particle(particle_current.begin() + begin, particle_current.begin() + end);
      Correlations block_corr =
          compute_correlations(block_heat, block_particle, config.ensemble_size, config.max_correlation_lag);
      RunningKappa block_running =
          integrate_correlations(block_corr, sample_dt, config.temperature, config.carrier_density_2d);
      if (!block_running.kappa.empty()) {
        block_kappa.push_back(block_running.kappa.back());
        block_kappa_j0.push_back(block_running.kappa_j0.back());
      }
    }
  }
  result.kappa_std_2d = stddev(block_kappa);
  result.kappa_j0_std_2d = stddev(block_kappa_j0);

  std::filesystem::create_directories(config.output_dir);
  std::string suffix = material.thermal_axis == 'y' ? "y" : "x";
  std::ofstream corr_file(config.output_dir + "/gk_heat_current_correlation_" + suffix + ".txt");
  corr_file << "# tau_s C_QQ C_QJ C_JJ\n";
  for (size_t lag = 0; lag < corr.qq.size(); ++lag) {
    corr_file << lag * sample_dt / units::s << " "
              << corr.qq[lag] << " "
              << corr.qj[lag] << " "
              << corr.jj[lag] << "\n";
  }

  std::ofstream kappa_file(config.output_dir + "/gk_kappa_running_" + suffix + ".txt");
  kappa_file << "# tau_s I_QQ I_QJ I_JJ kappa_W_per_K kappa_J0_W_per_K\n";
  for (size_t lag = 0; lag < running.kappa.size(); ++lag) {
    kappa_file << lag * sample_dt / units::s << " "
               << running.iqq[lag] << " "
               << running.iqj[lag] << " "
               << running.ijj[lag] << " "
               << running.kappa[lag] << " "
               << running.kappa_j0[lag] << "\n";
  }

  std::ofstream block_file(config.output_dir + "/gk_kappa_blocks_" + suffix + ".txt");
  block_file << "# block kappa_W_per_K kappa_J0_W_per_K\n";
  for (size_t i = 0; i < block_kappa.size(); ++i) {
    block_file << i << " " << block_kappa[i] << " " << block_kappa_j0[i] << "\n";
  }

  return result;
}

} // namespace probable
