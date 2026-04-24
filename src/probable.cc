#include <iostream>
#include <random>
#if __has_include(<omp.h>)
#include <omp.h>
#endif

#include <probable.hh>

namespace probable {

Vec2 Material::create_particle() const {
  double T_avg = (T_left + T_right) / 2;
  double m_eff = (mx + my) / 2;
  double p_max = 5 * sqrt(2 * m_eff * consts::kB * T_avg);
  while (true) {
    double prob = uniform() * exp(-Delta / (consts::kB * T_avg));
    Vec2 p = {p_max * (-1 + uniform() * 2), p_max * (-1 + uniform() * 2)};
    if (prob < exp(-energy(p) / (consts::kB * T_avg))) {      
      return p;
    }
  }
}

Pos2D Material::create_initial_position() const {
  return {uniform() * Lx, uniform() * Ly};  // равномерно от 0 до Lx
}

void Material::apply_boundary(double &x, Vec2 &p) const {
  double T = -1;
  double sign = 0;
  
  if (x < 0) {
    x = 0;
    T = T_left;
    sign = 1;
  } else if (x > Lx) {
    x = Lx;
    T = T_right;
    sign = -1;
  }
    
  if (T != -1) {
    double m_eff = (mx + my) / 2;
    double p_max = 5 * sqrt(2 * m_eff * consts::kB * T);
    double v_max = std::sqrt(Delta / mx);
    double boltzmann_floor = std::exp(-Delta / (consts::kB * T));
        
    while (true) {
      double p_x_new = p_max * sign * uniform();
      double p_y_new = p_max * (-1 + uniform() * 2);
      Vec2 p_test{p_x_new, p_y_new};
      double incoming_vx = sign * velocity(p_test).x;
      
      if (incoming_vx <= 0) {
        continue;
      }
      double prob = uniform() * v_max * boltzmann_floor;
      double weight = incoming_vx * std::exp(-energy(p_test) / (consts::kB * T));
      if (prob < weight) {
        p = p_test;
        break;
      }
    }
  }
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
    if (flags & DumpFlags::heat_flux) {
      heat_flux.push_back(q_flux);
    }
    if (flags & DumpFlags::particle_flux) {
      particle_flux.push_back(n_flux);
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
    if (r.flags & DumpFlags::heat_flux) {
      s << r.heat_flux[i] << " ";
    }
    if (r.flags & DumpFlags::particle_flux) {
      s << r.particle_flux[i] << " ";
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
#pragma omp parallel for
  for (size_t i = 0; i < ensemble_size; ++i) {
    Results &result = results[i];
    if (flags != DumpFlags::none) {
      result = Results(alloc, flags);
    }
    result.average_velocity = {0, 0};
    result.scattering_count.assign(mechanisms.size(), 0);
    Pos2D r = material.create_initial_position();
    Vec2 p = material.create_particle();
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

      r.y += v.y * time_step;
      r.x += v.x * time_step;
      if (r.y < 0) r.y += material.Ly;
      if (r.y >= material.Ly) r.y -= material.Ly;
      material.apply_boundary(r.x, p);
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
      result.append(j, j * time_step, p_, v, e, scattering_mechanism, x_old, heat, v.x);
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
