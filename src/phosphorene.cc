#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <probable.hh>

using namespace probable;
const double density = 1.4e-3 * units::g / pow(units::m, 2);     // density
const double sound_velocity = 6.8e3 * units::m / units::s;      // sound velocity
const double acoustic_deformation_potential = 4.9 * units::eV; // acoustic deformation potential, средний по двум напрпавлениям
const double Delta = 1.0 * units::eV; // полуширина запрещённой зоны
const double Lx = 1.0 * units::um; // размер области
const int Nx = 100 // количество бинов

struct AcousticScattering : public Scattering {
  double constant;  // (8π D_a² kT)/(h² ℏ ρ s²)
  double Da, s;  // параметры материала
  
  AcousticScattering(const Material &m, double temperature,
                       double Da_, double s_)
      : Scattering(m, 0), Da(Da_), s(s_) {
    double kT = consts::kB * temperature;
    
    // W_ac = (2 D_a² kT E)/(π ℏ³ ρ s²) * √((2 m_x m_y)/Δ(Δ+E)) * K((Δ-E)/(Δ+E))
    // Константа без E (E будет в rate())
    constant = 2 * Da * Da * kT / (math::pi * pow(consts::hbar, 3) * density * s * s);
  }
    
  double rate(const Vec2 &p) const override {
    double E = m.energy(p);  // энергия частицы (с учетом анизотропии и щели)
    
    double Dplus = m.Delta + E;
    double Dminus = m.Delta - E;
    
    // √((2 m_x m_y)/Δ(Δ+E))
    double sqrt_term = std::sqrt((2 * m.mx * m.my) / (m.Delta * Dplus));
    
    // Параметр для эллиптического интеграла k = (Δ-E)/(Δ+E)
    double k_param = Dminus / Dplus;
    
    // Полный эллиптический интеграл первого рода K(k)
    double K = std::comp_ellint_1(k_param);
    
    return constant * E * sqrt_term * K;
  }
};

// Испускание оптического фонона
struct OpticalEmissionScattering : public Scattering {
  double constant;
  double Do;
  
  OpticalEmissionScattering(const Material &m, double temperature,
                            double Do_, double phonon_energy)
      : Scattering(m, phonon_energy), Do(Do_) {
    
    double omega0 = phonon_energy / consts::hbar;
        
    // W_op = (D_o² E) / (π ℏ² ρ ω₀) * √((2 m_x m_y)/Δ(Δ+E)) * K((Δ-E)/(Δ+E))
    // Константа без E (E будет в rate())
    constant = Do * Do / (math::pi * std::pow(consts::hbar, 2) * density * omega0);
  }
  
  double rate(const Vec2 &p) const override {
    double E = m.energy(p) - energy;
    
    if (E < m.Delta) return 0;
    
    double Dplus = m.Delta + E;
    double Dminus = m.Delta - E;
    
    double sqrt_term = std::sqrt((2 * m.mx * m.my) / (m.Delta * Dplus));
    double k_param = Dminus / Dplus;
    double K = std::comp_ellint_1(k_param);
    
    return constant * E * sqrt_term * K;
  }
};

template <typename T> T parse(const std::string &s) {
  std::stringstream ss(s);
  T result;
  ss >> result;
  return result;
}

template <typename T> std::ostream &operator<<(std::ostream &s, std::vector<T> t) {
  s << "[";
  for (std::size_t i = 0; i < t.size(); i++) {
    s << t[i] << (i == t.size() - 1 ? "" : ", ");
  }
  return s << "]";
}

template <typename T> T sum(std::vector<T> t) {
  T s{};
  for (auto &v : t) {
    s += v;
  }
  return s;
}

int main(int argc, char const *argv[]) {
  if (argc != 7) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0] << " <ensemble size> <temperature> <Ex> <Ey> <Bz> <all_time>\n";
    return 1;
  }

  size_t ensemble_size = parse<int>(argv[1]);
  double temperature = parse<double>(argv[2]) * units::K;
  Vec2 electric_field{parse<double>(argv[3]) * units::V / units::m, 
                      parse<double>(argv[4]) * units::V / units::m};
  double magnetic_field_z = parse<double>(argv[5]) * units::T;
  double all_time = parse<double>(argv[6]) * units::s;

  Material phosphorene{
    1.285 * consts::me, // mx (ZZ)
    0.125 * consts::me, // my (AC)
    Delta,
    Lx,
    Nx
  };
  
  // Вектор механизмов рассеяния
  std::vector<Scattering *> mechanisms{
    new AcousticScattering(
      phosphorene, 
      temperature,
      acoustic_deformation_potential,
      sound_velocity
    ),
    new OpticalEmissionScattering(
      phosphorene,
      temperature,
      5.67e10 * units::eV / units::m, // TO2
      56e-3 * units::eV
    )
  };

  double time_step = 1e-16 * units::s;

  // print info on stdout
  std::cout << "Ensemble size:    " << ensemble_size << "\n";
  std::cout << "Sample length:    " << Lx / units::m << "\n";
  std::cout << "Time step:        " << time_step / units::s << " s\n";
  std::cout << "Simulation time:  " << all_time / units::s << " s\n";
  std::cout << "Temperature:      " << temperature / units::K << " K\n";
  std::cout << "Electric field:   " << electric_field / units::V * units::m << " V/m\n";
  std::cout << "Magnetic field:   " << magnetic_field_z / units::T << " T\n";
  std::cout << "Scattering mechanisms:\n";
  for (size_t i = 0; i < mechanisms.size(); ++i) {
    std::cout << "  " << i + 1 << ": " << *mechanisms[i] << '\n';
  }
  auto results = simulate(phosphorene,
                          mechanisms,
                          temperature,
                          electric_field,
                          magnetic_field_z,
                          time_step,
                          all_time,
                          ensemble_size,
                          DumpFlags(DumpFlags::none));
  // Обработка результатов
  Vec2 average_velocity;
  Vec2 average_velocity2;
  std::vector<double> scattering_rates(mechanisms.size(), 0);
  std::vector<double> total_counts(mechanisms.size(), 0);
  
  for (std::size_t i = 0; i < results.size(); ++i) {
    average_velocity += (results[i].average_velocity - average_velocity) / (i + 1);
    average_velocity2 +=
        (results[i].average_velocity * results[i].average_velocity - average_velocity2) / (i + 1);
        
    for (std::size_t j = 0; j < mechanisms.size(); ++j) {
      scattering_rates[j] +=
          (results[i].scattering_count[j] / all_time * units::s - scattering_rates[j]) / (i + 1);
	  total_counts[j] += results[i].scattering_count[j];
    }
  }
  Vec2 std_velocity = (average_velocity2 - average_velocity * average_velocity).sqrt();
  
  std::cout << "\n===== Results =====\n";
  std::cout << "      Directions: {ZZ, AC}\n";
  std::cout << "Average velocity: " << average_velocity << " μm/ps\n";
  std::cout << "             std: " << std_velocity << " μm/ps\n";
  std::cout << "Scattering rates: " << scattering_rates << " 1/s\n";
  std::cout << "           Total: " << sum(scattering_rates) << " 1/s\n";
  std::cout << "Scattering count: " << total_counts << " times\n";
  
  return 0;
}
