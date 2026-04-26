#pragma once

#include <cstdint>
#include <cxxabi.h>
#include <cmath>
#include <ostream>
#include <typeinfo>
#include <vector>

#include <vec2.hh>
#include <position.hh>

namespace units {
// energy
const double eV = 1;
const double J = 1 / 1.6e-19;

// length
const double um = 1;
const double cm = um * 1e4;
const double m = um * 1e6;
const double nm = um * 1e-3;

// time
const double ps = 1;
const double s = ps * 1e12;

// mass
const double kg = J * s * s / m / m;
const double g = 1e-3 * kg;

// potential
const double V = 1;

// magnetic field
const double T = 3e8 * (V / m) / (m / s);

// electric charge
const double C = J / V;

// temperature
const double K = 1.38e-23 / 1.6e-19;
} // namespace units

namespace consts {
const double hbar = 1.05e-34 * units::J * units::s;
const double e = 1.6e-19 * units::C;
const double c = 3e8 * units::m / units::s;
const double me = 9.1e-31 * units::kg;
const double kB = 1.38e-23 * units::J / units::K;
const double eps0 = 8.85e-12 * units::C / units::V / units::m;
} // namespace consts

namespace math {
const double e = exp(1);
const double pi = acos(-1);

inline double comp_ellint_1(double k) {
  // Complete elliptic integral K(k) via the arithmetic-geometric mean.
  // This keeps the code buildable with libc++ versions lacking std::comp_ellint_1.
  double kk = std::fabs(k);
  if (kk >= 1) {
    kk = std::nextafter(1.0, 0.0);
  }
  double a = 1.0;
  double b = std::sqrt(1.0 - kk * kk);
  for (int i = 0; i < 32 && std::fabs(a - b) > 1e-14 * a; ++i) {
    double next_a = 0.5 * (a + b);
    b = std::sqrt(a * b);
    a = next_a;
  }
  return pi / (2.0 * a);
}
} // namespace math

namespace probable {

struct Material {
  double mx, my;             // эффективные массы по осям
  double Delta;              // полуширина запрещённой зоны
  double Lx, Ly;             // размеры области
  int Nx;                    // количество ячеек
  double cell_size;          // Lx / Nx
  double T_left, T_right;    // Температура на краях
  
  Material(double mx_, double my_, double Delta_ = 0, double Lx_ = 0, 
           double Ly_ = 0, int Nx_ = 1, double T_l = 0, double T_r = 0) 
    : mx(mx_), my(my_), Delta(Delta_), Lx(Lx_), Ly(Ly_), Nx(Nx_), 
      cell_size(Lx_ / Nx_), T_left(T_l), T_right(T_r) {}
  
  double energy(const Vec2 &p) const {
    // ϵ(p) = √((Δ p_x²)/m_x + (Δ + p_y²/(2m_y))²)
    double term1 = Delta * p.x * p.x / mx;
    double term2 = Delta + p.y * p.y / (2 * my);
    return std::sqrt(term1 + term2 * term2);
  }
  
  Vec2 velocity(const Vec2 &p) const {
    double vx = (Delta * p.x) / (mx * std::sqrt(Delta / mx * p.x * p.x + std::pow((Delta + p.y * p.y / (2 * my)), 2)));
    double vy = p.y * (Delta + p.y * p.y / (2 * my)) / (my * std::sqrt(Delta / mx * p.x * p.x + std::pow((Delta + p.y * p.y / (2 * my)), 2)));
    return {vx, vy};
  }
  
  Vec2 create_particle() const;

  Vec2 create_particle_at_temperature(double T) const;

  Vec2 create_flux_particle_at_temperature(double T, double sign) const;
  
  Pos2D create_initial_position() const;
  
  // Отражающие граничные условия
  int apply_boundary(double &x, Vec2 &p) const;
  
  // получение индекса ячейки
  int get_cell_index(double x) const {
    if (x <= 0) return 0;
    if (x >= Lx) return Nx - 1;
    return static_cast<int>(x / cell_size);
  }
  
  double get_temperature(double x) const {
    return T_left + (T_right - T_left) * (x / Lx);
  }

  double mean_excess_energy_at_temperature(double T) const;

  double mean_flux_excess_energy_at_temperature(double T) const;
  
  double temperature_from_mean_excess_energy(double mean_excess, double T_min, double T_max) const;

  double temperature_from_mean_flux_excess_energy(double mean_excess, double T_min, double T_max) const;
};

struct Scattering {
  const Material &m;
  const double energy;
  Scattering(const Material &m, double e) : m(m), energy(e) {}
  virtual double rate(const Vec2 &p, double x) const = 0;
  Vec2 scatter(const Vec2 &p) const;
  virtual ~Scattering() {}
};

inline std::ostream &operator<<(std::ostream &s, const Scattering &sc) {
  int status;
  char *realname = abi::__cxa_demangle(typeid(sc).name(), 0, 0, &status);
  s << realname << " " << sc.energy;
  free(realname);
  return s;
}

enum DumpFlags {
  // contents
  none = 0,
  number = 1,
  time = number << 1,
  momentum = time << 1,
  energy = momentum << 1,
  position = energy << 1,
  velocity = position << 1,
  scattering = velocity << 1,
  heat_flux = scattering << 1,
  particle_flux = heat_flux << 1,
  all = number | time | momentum | energy | velocity | position | scattering | heat_flux | particle_flux,
  // frequency
  // without this flag it will dump on every step
  on_scatterings = particle_flux << 1,
};

struct Results {
  size_t size;
  DumpFlags flags;
  Vec2 average_velocity;
  std::vector<uint32_t> scattering_count;
  std::vector<uint32_t> ns;
  std::vector<double> ts;
  std::vector<Vec2> momentums;
  std::vector<Vec2> velocities;
  std::vector<double> energies;
  std::vector<uint32_t> scatterings;
  std::vector<double> positions;
  double heat_flux_sum;
  double particle_flux_sum;
  uint64_t flux_samples;
  size_t flux_sample_stride;
  std::vector<double> heat_flux_windows;
  std::vector<double> particle_flux_windows;
  std::vector<uint64_t> flux_window_samples;
  std::vector<double> bin_excess_energy;
  std::vector<uint64_t> bin_samples;
  std::vector<double> initial_bin_excess_energy;
  std::vector<uint64_t> initial_bin_samples;
  double injected_left_excess_energy;
  double injected_right_excess_energy;
  uint64_t injected_left_samples;
  uint64_t injected_right_samples;
  Results() {}
  Results(size_t cap, DumpFlags flags = DumpFlags::none, size_t flux_windows = 0, size_t flux_stride = 1)
      : size(0), flags(flags), heat_flux_sum(0), particle_flux_sum(0), flux_samples(0),
        flux_sample_stride(flux_stride), injected_left_excess_energy(0), injected_right_excess_energy(0),
        injected_left_samples(0), injected_right_samples(0),
        ns(), ts(), momentums(), velocities(), energies(), scatterings() {
    ns.reserve(cap);
    ts.reserve(cap);
    momentums.reserve(cap);
    velocities.reserve(cap);
    energies.reserve(cap);
    scatterings.reserve(cap);
    positions.reserve(cap);
    heat_flux_windows.assign(flux_windows, 0);
    particle_flux_windows.assign(flux_windows, 0);
    flux_window_samples.assign(flux_windows, 0);
  }
  void append(uint32_t n, double t, const Vec2 &p, const Vec2 &v, double e, size_t s, double x, double q_flux, double n_flux);
  friend std::ostream &operator<<(std::ostream &s, const Results &r);
};

std::vector<Results> simulate(const Material &material,
                              const std::vector<Scattering *> mechanisms,
                              const Vec2 &electric_field,
                              double magnetic_field,
                              double time_step,
                              double all_time,
                              size_t ensemble_size,
                              DumpFlags flags = DumpFlags::all);
                              
// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable
