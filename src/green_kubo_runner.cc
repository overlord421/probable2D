#include <green_kubo_runner.hh>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
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
            p = mechanisms[k]->scatter(p);
            free_flight[k] = -std::log(uniform());
            local_scattering_counts[k] += 1;
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
