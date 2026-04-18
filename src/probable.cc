#include <iostream>
#include <omp.h>
#include <random>

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

Pos1D Material::create_initial_position() const {
  return {uniform() * Lx};  // равномерно от 0 до Lx
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
    
    double E = energy(p);
    double Dplus = Delta + E;
    double Dminus = Delta - E;
    double g = 2; // спиновое вырождение
    double sqrt_term = std::sqrt((2 * mx * my) / (Delta * Dplus));
    // Параметр для эллиптического интеграла k = (Δ-E)/(Δ+E)
    double k_param = Dminus / Dplus;
    // Полный эллиптический интеграл первого рода K(k)
    double K = std::comp_ellint_1(k_param);
    // плотность состояний (DOS)
    double rho = g * E * sqrt_term * K / pow(math::pi * consts::hbar, 2);
        
    while (true) {
      double p_x_new = p_max * sign * uniform();
      double prob = uniform() * exp(-Delta / (consts::kB * T));
      
      Vec2 p_test{p_x_new, p.y};
      if (prob < exp(-energy(p_test) / (consts::kB * T)) * fabs(p.x) * rho) {
        p.x = p_x_new;
        break;
      }
    }
  }
}

Vec2 Scattering::scatter(const Vec2 &p) const {
  // В 2D только угол φ
  double phi = 2 * math::pi * uniform();
  
  double e = m.energy(p) - energy;
  // нахожу p0, решая квадратное уравнение
  double a = pow(sin(phi), 4) / (4 * pow(m.my, 2));
  double b = m.Delta * (pow(cos(phi), 2) / m.mx + pow(sin(phi), 2) / m.my);
  double c = -(pow(e, 2) - pow(m.Delta, 2));
  double D = pow(b, 2) - (4 * a * c);
  double p0 = sqrt((sqrt(D) - b) / (2 * a));
  
  return {p0 * cos(phi), p0 * sin(phi)};
}

void Results::append(uint32_t n, double t, const Vec2 &p, const Vec2 &v, double e, size_t s, double x, double flux) {
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
    if (flags & DumpFlags::energy_flux) {
      energy_flux.push_back(flux);
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
    if (r.flags & DumpFlags::energy_flux) {
      s << r.energy_flux[i] << " ";
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
    Pos1D x = material.create_initial_position();
    Vec2 p = material.create_particle();
    std::vector<double> free_flight(mechanisms.size(), 0);
    for (double &l : free_flight) {
      l = -log(uniform());
    }

    for (size_t j = 0; j < steps; ++j) {
      Vec2 p_ = p;
      Vec2 v = material.velocity(p_);
      double e = material.energy(p_);
      size_t scattering_mechanism = 0; // means no scattering
      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p_, x.x) * time_step;
        if (free_flight[k] < 0) {
          p = mechanisms[k]->scatter(p_);
          free_flight[k] = -log(uniform());
          scattering_mechanism = k + 1; // enumerate mechanisms from 1
          break;
        }
      }
      
      double flux = e * v.x;
      result.append(j, j * time_step, p_, v, e, scattering_mechanism, x.x, flux);      
      if (not scattering_mechanism) {
        x.x += v.x * time_step;  // обновляем позицию
        material.apply_boundary(x.x, p);  // применяем граничные условия
        // Для B, направленного перпендикулярно плоскости (только Bz):
        p += -consts::e * (electric_field + v.cross_with_B(magnetic_field_z)) * time_step;
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