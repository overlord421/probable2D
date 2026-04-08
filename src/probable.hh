#pragma once

#include <cstdint>
#include <cxxabi.h>
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
} // namespace math

namespace probable {

struct Material {
  double mx, my;        // эффективные массы по осям
  double Delta;         // полуширина запрещённой зоны
  double Lx;            // размер области по x
  int Nx;               // количество ячеек
  double cell_size;     // Lx / Nx
  
  Material(double mx_, double my_, double Delta_ = 0, double Lx_ = 0, int Nx_ = 1) 
    : mx(mx_), my(my_), Delta(Delta_), Lx(Lx_), Nx(Nx_), cell_size(Lx_ / Nx_) {}
  
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
  
  Vec2 create_particle(double temperature) const;
  
  Pos1D create_initial_position() const {
    return {uniform() * Lx};  // равномерно от 0 до Lx
  }
  
  // Отражающие граничные условия
  double apply_boundary(double x) const {
    if (x < 0) return -x;           // отражение от левой границы
    if (x > Lx) return 2*Lx - x;    // отражение от правой границы
    return x;
  }
  
  // получение индекса ячейки
  int get_cell_index(double x) const {
    if (x <= 0) return 0;
    if (x >= Lx) return Nx - 1;
    return static_cast<int>(x / cell_size);
  }
};

struct Scattering {
  const Material &m;
  const double energy;
  Scattering(const Material &m, double e) : m(m), energy(e) {}
  virtual double rate(const Vec2 &p) const = 0;
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
  all = number | time | momentum | energy | velocity | position | scattering,
  // frequency
  // without this flag it will dump on every step
  on_scatterings = scattering << 1,
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
  Results() {}
  Results(size_t cap, DumpFlags flags = DumpFlags::none)
      : size(0), flags(flags), ns(), ts(), momentums(), velocities(), energies(), scatterings() {
    ns.reserve(cap);
    ts.reserve(cap);
    momentums.reserve(cap);
    velocities.reserve(cap);
    energies.reserve(cap);
    scatterings.reserve(cap);
  }
  void append(uint32_t n, double t, const Vec2 &p, const Vec2 &v, double e, size_t s, double x);
  friend std::ostream &operator<<(std::ostream &s, const Results &r);
};

std::vector<Results> simulate(const Material &material,
                              const std::vector<Scattering *> mechanisms,
                              double temperature,
                              const Vec2 &electric_field,
                              double magnetic_field,
                              double time_step,
                              double all_time,
                              size_t ensemble_size,
                              DumpFlags flags = DumpFlags::all);
                              
// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable