#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#if __has_include(<omp.h>)
#include <omp.h>
#endif

#include <probable.hh>

namespace probable {

namespace {

double fermi_occupation_at(double E, double mu, double T) {
  double x = (E - mu) / (consts::kB * T);
  if (x > 50) {
    return 0;
  }
  if (x < -50) {
    return 1;
  }
  return 1 / (std::exp(x) + 1);
}

double fermi_density_at_mu(const Material &material, double T, double mu) {
  double m_eff = (material.mx + material.my) / 2;
  double excess_scale = std::max(consts::kB * T, mu - material.Delta + 16 * consts::kB * T);
  double p_max = 6 * std::sqrt(2 * m_eff * excess_scale);
  const int grid = 96;
  double dp = 2 * p_max / grid;
  double occupation_sum = 0;

  for (int ix = 0; ix < grid; ++ix) {
    double px = p_max * (-1 + 2.0 * (ix + 0.5) / grid);
    for (int iy = 0; iy < grid; ++iy) {
      double py = p_max * (-1 + 2.0 * (iy + 0.5) / grid);
      occupation_sum += fermi_occupation_at(material.energy({px, py}), mu, T);
    }
  }

  double state_density = material.spin_degeneracy / std::pow(2 * math::pi * consts::hbar, 2);
  return state_density * occupation_sum * dp * dp;
}

} // namespace

Vec2 Material::create_particle() const {
  double T_avg = (T_left + T_right) / 2;
  return create_particle_at_temperature(T_avg);
}

Vec2 Material::create_particle_at_temperature(double T) const {
  double m_eff = (mx + my) / 2;
  double mu = chemical_potential(T);
  double excess_scale = use_fermi_dirac ? std::max(consts::kB * T, mu - Delta + 8 * consts::kB * T)
                                        : consts::kB * T;
  double p_max = 5 * sqrt(2 * m_eff * excess_scale);
  while (true) {
    Vec2 p = {p_max * (-1 + uniform() * 2), p_max * (-1 + uniform() * 2)};
    double accept = use_fermi_dirac ? occupation(energy(p), T)
                                    : std::exp(-(energy(p) - Delta) / (consts::kB * T));
    if (uniform() < accept) {
      return p;
    }
  }
}

Vec2 Material::create_flux_particle_at_temperature(double T, double sign) const {
  double m_eff = (mx + my) / 2;
  double mu = chemical_potential(T);
  double excess_scale = use_fermi_dirac ? std::max(consts::kB * T, mu - Delta + 8 * consts::kB * T)
                                        : consts::kB * T;
  double p_max = 5 * sqrt(2 * m_eff * excess_scale);
  double v_max = 0;
  for (int i = 1; i <= 64; ++i) {
    Vec2 p_probe;
    if (thermal_axis == 'y') {
      p_probe = {0, p_max * sign * i / 64.0};
    } else {
      p_probe = {p_max * sign * i / 64.0, 0};
    }
    v_max = std::max(v_max, sign * velocity_component(velocity(p_probe)));
  }

  while (true) {
    Vec2 p;
    if (thermal_axis == 'y') {
      p = {p_max * (-1 + uniform() * 2), p_max * sign * uniform()};
    } else {
      p = {p_max * sign * uniform(), p_max * (-1 + uniform() * 2)};
    }
    double incoming_v = sign * velocity_component(velocity(p));
    if (incoming_v <= 0) {
      continue;
    }

    double occupation_weight = use_fermi_dirac ? occupation(energy(p), T)
                                               : std::exp(-(energy(p) - Delta) / (consts::kB * T));
    double accept = incoming_v / v_max * occupation_weight;
    if (uniform() < accept) {
      return p;
    }
  }
}

Pos2D Material::create_initial_position() const {
  return {uniform() * Lx, uniform() * Ly};  // равномерно от 0 до Lx
}

int Material::apply_boundary(Pos2D &r, Vec2 &p) const {
  double T = -1;
  double sign = 0;
  int side = 0;
  double coord = axis_coordinate(r);
  
  if (coord < 0) {
    T = T_left;
    sign = 1;
    side = -1;
    if (thermal_axis == 'y') {
      r.y = 0;
    } else {
      r.x = 0;
    }
  } else if (coord > axis_length()) {
    T = T_right;
    sign = -1;
    side = 1;
    if (thermal_axis == 'y') {
      r.y = Ly;
    } else {
      r.x = Lx;
    }
  }
    
  if (T != -1) {
    p = create_flux_particle_at_temperature(T, sign);
  }
  return side;
}

double Material::effective_density_of_states_2d(double T) const {
  double density_mass = std::sqrt(mx * my);
  return spin_degeneracy * density_mass * consts::kB * T / (2 * math::pi * consts::hbar * consts::hbar);
}

double Material::chemical_potential(double T) const {
  if (!use_fermi_dirac || carrier_density_2d <= 0) {
    return -1e6 * units::eV;
  }

  long long temperature_key = std::llround(T / units::K * 100);
  using CacheKey = std::pair<const Material *, long long>;
  thread_local std::map<CacheKey, double> cache;
  CacheKey key{this, temperature_key};
  auto found = cache.find(key);
  if (found != cache.end()) {
    return found->second;
  }

  double Nc = effective_density_of_states_2d(T);
  double eta_guess = carrier_density_2d / Nc > 50
      ? carrier_density_2d / Nc
      : std::log(std::expm1(carrier_density_2d / Nc));
  double mu_guess = Delta + consts::kB * T * eta_guess;
  double low = Delta - 80 * consts::kB * T;
  double high = std::max(Delta + 80 * consts::kB * T, mu_guess + 40 * consts::kB * T);

  for (int i = 0; i < 16 && fermi_density_at_mu(*this, T, high) < carrier_density_2d; ++i) {
    high += (high - low);
  }

  for (int i = 0; i < 48; ++i) {
    double mid = 0.5 * (low + high);
    if (fermi_density_at_mu(*this, T, mid) < carrier_density_2d) {
      low = mid;
    } else {
      high = mid;
    }
  }

  double mu = 0.5 * (low + high);
  cache[key] = mu;
  return mu;
}

double Material::occupation(double E, double T) const {
  if (!use_fermi_dirac) {
    return 0;
  }
  return fermi_occupation_at(E, chemical_potential(T), T);
}

double Material::mean_excess_energy_at_temperature(double T) const {
  double m_eff = (mx + my) / 2;
  double mu = chemical_potential(T);
  double excess_scale = use_fermi_dirac ? std::max(consts::kB * T, mu - Delta + 8 * consts::kB * T)
                                        : consts::kB * T;
  double p_max = 6 * sqrt(2 * m_eff * excess_scale);
  const int grid = 72;
  double weighted_energy = 0;
  double weight_sum = 0;

  for (int ix = 0; ix < grid; ++ix) {
    double px = p_max * (-1 + 2.0 * (ix + 0.5) / grid);
    for (int iy = 0; iy < grid; ++iy) {
      double py = p_max * (-1 + 2.0 * (iy + 0.5) / grid);
      Vec2 p{px, py};
      double E = energy(p);
      double weight = use_fermi_dirac ? occupation(E, T) : std::exp(-E / (consts::kB * T));
      weighted_energy += (E - Delta) * weight;
      weight_sum += weight;
    }
  }

  if (weight_sum == 0) {
    return 0;
  }
  return weighted_energy / weight_sum;
}

double Material::mean_flux_excess_energy_at_temperature(double T) const {
  double m_eff = (mx + my) / 2;
  double mu = chemical_potential(T);
  double excess_scale = use_fermi_dirac ? std::max(consts::kB * T, mu - Delta + 8 * consts::kB * T)
                                        : consts::kB * T;
  double p_max = 6 * sqrt(2 * m_eff * excess_scale);
  const int grid = 72;
  double weighted_energy = 0;
  double weight_sum = 0;

  for (int ix = 0; ix < grid; ++ix) {
    double px = thermal_axis == 'y'
        ? p_max * (-1 + 2.0 * (ix + 0.5) / grid)
        : p_max * (ix + 0.5) / grid;
    for (int iy = 0; iy < grid; ++iy) {
      double py = thermal_axis == 'y'
          ? p_max * (iy + 0.5) / grid
          : p_max * (-1 + 2.0 * (iy + 0.5) / grid);
      Vec2 p{px, py};
      double E = energy(p);
      double occupation_weight = use_fermi_dirac ? occupation(E, T) : std::exp(-E / (consts::kB * T));
      double weight = velocity_component(velocity(p)) * occupation_weight;
      weighted_energy += (E - Delta) * weight;
      weight_sum += weight;
    }
  }

  if (weight_sum == 0) {
    return 0;
  }
  return weighted_energy / weight_sum;
}

double Material::temperature_from_mean_excess_energy(double mean_excess, double T_min, double T_max) const {
  if (mean_excess <= 0) {
    return T_min;
  }

  double high = T_max;
  for (int i = 0; i < 8 && mean_excess_energy_at_temperature(high) < mean_excess; ++i) {
    high *= 2;
  }

  double low = T_min;
  for (int i = 0; i < 48; ++i) {
    double mid = 0.5 * (low + high);
    if (mean_excess_energy_at_temperature(mid) < mean_excess) {
      low = mid;
    } else {
      high = mid;
    }
  }
  return 0.5 * (low + high);
}

double Material::temperature_from_mean_flux_excess_energy(double mean_excess, double T_min, double T_max) const {
  if (mean_excess <= 0) {
    return T_min;
  }

  double high = T_max;
  for (int i = 0; i < 8 && mean_flux_excess_energy_at_temperature(high) < mean_excess; ++i) {
    high *= 2;
  }

  double low = T_min;
  for (int i = 0; i < 48; ++i) {
    double mid = 0.5 * (low + high);
    if (mean_flux_excess_energy_at_temperature(mid) < mean_excess) {
      low = mid;
    } else {
      high = mid;
    }
  }
  return 0.5 * (low + high);
}

Vec2 Scattering::scatter(const Vec2 &p) const {
  // В 2D только угол φ
  double phi = 2 * math::pi * uniform();
  
  double e = m.energy(p) - energy;
  if (e <= m.Delta) {
    return p;
  }
  // нахожу p0, решая квадратное уравнение
  double a = pow(sin(phi), 4) / (4 * pow(m.my, 2));
  double b = m.Delta * (pow(cos(phi), 2) / m.mx + pow(sin(phi), 2) / m.my);
  double c = -(pow(e, 2) - pow(m.Delta, 2));
  if (std::fabs(a) < 1e-300) {
    double p0 = std::sqrt(-c / b);
    return {p0 * cos(phi), p0 * sin(phi)};
  }
  double D = pow(b, 2) - (4 * a * c);
  double p0 = sqrt((sqrt(D) - b) / (2 * a));
  
  return {p0 * cos(phi), p0 * sin(phi)};
}

void Results::append(uint32_t n, double t, const Vec2 &p, const Vec2 &v, double e, size_t s, double x, double q_flux, double n_flux) {
  average_velocity += (v - average_velocity) / (n + 1);
  if (s) {
    scattering_count[s - 1] += 1;
  }
  if (flags & DumpFlags::heat_flux) {
    heat_flux_sum += q_flux;
  }
  if (flags & DumpFlags::particle_flux) {
    particle_flux_sum += n_flux;
  }
  if (flags & (DumpFlags::heat_flux | DumpFlags::particle_flux)) {
    flux_samples += 1;
    size_t window = flux_sample_stride > 0 ? n / flux_sample_stride : 0;
    if (window < flux_window_samples.size()) {
      heat_flux_windows[window] += q_flux;
      particle_flux_windows[window] += n_flux;
      flux_window_samples[window] += 1;
    }
  }
  if (not(flags & DumpFlags::on_scatterings) or s) {
    if (flags & DumpFlags::scattering) {
      scatterings.push_back(s);
    }
    if (flags & DumpFlags::number) {
      ns.push_back(n);
    }
    if (flags & DumpFlags::time) {
      ts.push_back(t);
    }
    if (flags & DumpFlags::momentum) {
      momentums.push_back(p);
    }
    if (flags & DumpFlags::velocity) {
      velocities.push_back(v);
    }
    if (flags & DumpFlags::energy) {
      energies.push_back(e);
    }
    if (flags & DumpFlags::position) {
      positions.push_back(x);
    }
    size += 1;
  }
}

std::ostream &operator<<(std::ostream &s, const Results &r) {
  if (r.flags == DumpFlags::none) {
    return s;
  }
  for (size_t i = 0; i < r.size; ++i) {
    if (r.flags & DumpFlags::number) {
      s << r.ns[i] << " ";
    }
    if (r.flags & DumpFlags::time) {
      s << r.ts[i] / units::s << " ";
    }
    if (r.flags & DumpFlags::momentum) {
      s << r.momentums[i].x << " ";
      s << r.momentums[i].y << " ";
      // Убрали .z компонент
    }
    if (r.flags & DumpFlags::velocity) {
      s << r.velocities[i].x / units::m * units::s << " ";
      s << r.velocities[i].y / units::m * units::s << " ";
      // Убрали .z компонент
    }
    if (r.flags & DumpFlags::energy) {
      s << r.energies[i] / units::eV << " ";
    }
    if (r.flags & DumpFlags::position) {
      s << r.positions[i] / units::m << " ";
    }
    if (r.flags & DumpFlags::scattering) {
      s << r.scatterings[i];
    }
    s << "\n";
  }
  return s;
}

std::vector<Results> simulate(const Material &material,
                              const std::vector<Scattering *> mechanisms,
                              const Vec2 &electric_field,
                              double magnetic_field_z,
                              double time_step,
                              double all_time,
                              size_t ensemble_size,
                              DumpFlags flags) {
  std::vector<Results> results(ensemble_size);
  size_t steps = all_time / time_step + 1;
  size_t alloc = (flags & DumpFlags::on_scatterings) ? steps / 10 : steps;
  if (flags == DumpFlags(DumpFlags::heat_flux | DumpFlags::particle_flux)) {
    alloc = 0;
  }
  size_t max_flux_windows = 1000;
  size_t adaptive_flux_stride = (steps + max_flux_windows - 1) / max_flux_windows;
  size_t flux_sample_stride = adaptive_flux_stride > 100 ? adaptive_flux_stride : 100;
  size_t flux_windows = (steps + flux_sample_stride - 1) / flux_sample_stride;
  size_t tally_start_step = steps / 5;
#pragma omp parallel for
  for (size_t i = 0; i < ensemble_size; ++i) {
    Results &result = results[i];
    if (flags != DumpFlags::none) {
      result = Results(alloc, flags, flux_windows, flux_sample_stride);
    }
    result.average_velocity = {0, 0};
    result.scattering_count.assign(mechanisms.size(), 0);
    result.bin_excess_energy.assign(material.Nx, 0);
    result.bin_samples.assign(material.Nx, 0);
    result.initial_bin_excess_energy.assign(material.Nx, 0);
    result.initial_bin_samples.assign(material.Nx, 0);
    Pos2D r = material.create_initial_position();
    Vec2 p = material.create_particle_at_temperature(material.get_temperature(material.axis_coordinate(r)));
    int initial_bin = material.get_cell_index(material.axis_coordinate(r));
    result.initial_bin_excess_energy[initial_bin] += material.energy(p) - material.Delta;
    result.initial_bin_samples[initial_bin] += 1;
    std::vector<double> free_flight(mechanisms.size(), 0);
    for (double &l : free_flight) {
      l = -log(uniform());
    }

    for (size_t j = 0; j < steps; ++j) {
      Vec2 p_ = p;
      Vec2 v = material.velocity(p_);
      double e = material.energy(p_);
      double coord_old = material.axis_coordinate(r);
      size_t scattering_mechanism = 0; // means no scattering
      double v_axis = material.velocity_component(v);
      double heat = (e - material.Delta) * v_axis;
      bool tally = j >= tally_start_step;
      if (tally) {
        int bin = material.get_cell_index(coord_old);
        result.bin_excess_energy[bin] += e - material.Delta;
        result.bin_samples[bin] += 1;
      }

      r.y += v.y * time_step;
      r.x += v.x * time_step;
      if (material.thermal_axis == 'y') {
        if (r.x < 0) r.x += material.Lx;
        if (r.x >= material.Lx) r.x -= material.Lx;
      } else {
        if (r.y < 0) r.y += material.Ly;
        if (r.y >= material.Ly) r.y -= material.Ly;
      }
      int injected_side = material.apply_boundary(r, p);
      if (injected_side < 0) {
        result.injected_left_excess_energy += material.energy(p) - material.Delta;
        result.injected_left_samples += 1;
      } else if (injected_side > 0) {
        result.injected_right_excess_energy += material.energy(p) - material.Delta;
        result.injected_right_samples += 1;
      }
      p += -consts::e * (electric_field + v.cross_with_B(magnetic_field_z)) * time_step;

      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p, material.axis_coordinate(r)) * time_step;
        if (free_flight[k] < 0) {
          Vec2 scattered = mechanisms[k]->scatter(p);
          double final_occupation =
              material.use_fermi_dirac ? material.occupation(material.energy(scattered), material.get_temperature(material.axis_coordinate(r))) : 0;
          bool accepted = !material.use_fermi_dirac || uniform() >= final_occupation;
          if (accepted) {
            p = scattered;
            scattering_mechanism = k + 1; // enumerate mechanisms from 1
          }
          free_flight[k] = -log(uniform());
          break;
        }
      }
      if (tally) {
        result.append(j, j * time_step, p_, v, e, scattering_mechanism, coord_old, heat, v_axis);
      } else {
        if (scattering_mechanism) {
          result.scattering_count[scattering_mechanism - 1] += 1;
        }
      }
    }
  }
  return results;
}

double uniform() {
  thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
}
} // namespace probable
