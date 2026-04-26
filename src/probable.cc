#include <iostream>
#include <random>
#if __has_include(<omp.h>)
#include <omp.h>
#endif

#include <probable.hh>

namespace probable {

Vec2 Material::create_particle() const {
  double T_avg = (T_left + T_right) / 2;
  return create_particle_at_temperature(T_avg);
}

Vec2 Material::create_particle_at_temperature(double T) const {
  double m_eff = (mx + my) / 2;
  double p_max = 5 * sqrt(2 * m_eff * consts::kB * T);
  while (true) {
    double prob = uniform() * exp(-Delta / (consts::kB * T));
    Vec2 p = {p_max * (-1 + uniform() * 2), p_max * (-1 + uniform() * 2)};
    if (prob < exp(-energy(p) / (consts::kB * T))) {      
      return p;
    }
  }
}

Vec2 Material::create_flux_particle_at_temperature(double T, double sign) const {
  double m_eff = (mx + my) / 2;
  double p_max = 5 * sqrt(2 * m_eff * consts::kB * T);
  double v_max = std::sqrt(Delta / mx);
  double boltzmann_floor = std::exp(-Delta / (consts::kB * T));

  while (true) {
    Vec2 p{p_max * sign * uniform(), p_max * (-1 + uniform() * 2)};
    double incoming_vx = sign * velocity(p).x;
    if (incoming_vx <= 0) {
      continue;
    }

    double prob = uniform() * v_max * boltzmann_floor;
    double weight = incoming_vx * std::exp(-energy(p) / (consts::kB * T));
    if (prob < weight) {
      return p;
    }
  }
}

Pos2D Material::create_initial_position() const {
  return {uniform() * Lx, uniform() * Ly};  // равномерно от 0 до Lx
}

int Material::apply_boundary(double &x, Vec2 &p) const {
  double T = -1;
  double sign = 0;
  int side = 0;
  
  if (x < 0) {
    x = 0;
    T = T_left;
    sign = 1;
    side = -1;
  } else if (x > Lx) {
    x = Lx;
    T = T_right;
    sign = -1;
    side = 1;
  }
    
  if (T != -1) {
    p = create_flux_particle_at_temperature(T, sign);
  }
  return side;
}

double Material::mean_excess_energy_at_temperature(double T) const {
  double m_eff = (mx + my) / 2;
  double p_max = 6 * sqrt(2 * m_eff * consts::kB * T);
  const int grid = 72;
  double weighted_energy = 0;
  double weight_sum = 0;

  for (int ix = 0; ix < grid; ++ix) {
    double px = p_max * (-1 + 2.0 * (ix + 0.5) / grid);
    for (int iy = 0; iy < grid; ++iy) {
      double py = p_max * (-1 + 2.0 * (iy + 0.5) / grid);
      Vec2 p{px, py};
      double E = energy(p);
      double weight = std::exp(-E / (consts::kB * T));
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
    Vec2 p = material.create_particle_at_temperature(material.get_temperature(r.x));
    int initial_bin = material.get_cell_index(r.x);
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
      double x_old = r.x;
      size_t scattering_mechanism = 0; // means no scattering
      double heat = (e - material.Delta) * v.x;
      bool tally = j >= tally_start_step;
      if (tally) {
        int bin = material.get_cell_index(x_old);
        result.bin_excess_energy[bin] += e - material.Delta;
        result.bin_samples[bin] += 1;
      }

      r.y += v.y * time_step;
      r.x += v.x * time_step;
      if (r.y < 0) r.y += material.Ly;
      if (r.y >= material.Ly) r.y -= material.Ly;
      int injected_side = material.apply_boundary(r.x, p);
      if (injected_side < 0) {
        result.injected_left_excess_energy += material.energy(p) - material.Delta;
        result.injected_left_samples += 1;
      } else if (injected_side > 0) {
        result.injected_right_excess_energy += material.energy(p) - material.Delta;
        result.injected_right_samples += 1;
      }
      p += -consts::e * (electric_field + v.cross_with_B(magnetic_field_z)) * time_step;

      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p, r.x) * time_step;
        if (free_flight[k] < 0) {
          p = mechanisms[k]->scatter(p);
          free_flight[k] = -log(uniform());
          scattering_mechanism = k + 1; // enumerate mechanisms from 1
          break;
        }
      }
      if (tally) {
        result.append(j, j * time_step, p_, v, e, scattering_mechanism, x_old, heat, v.x);
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
