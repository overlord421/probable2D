#include <fstream>
#include <filesystem>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <probable.hh>

using namespace probable;
const double density = 1.4e-3 * units::g / pow(units::m, 2);     // density
const double sound_velocity = 6.8e3 * units::m / units::s;      // sound velocity
const double acoustic_deformation_potential = 4.9 * units::eV; // acoustic deformation potential, средний по двум напрпавлениям
const double Delta = 1.0 * units::eV; // полуширина запрещённой зоны
const double Lx = 300 * units::nm; // размер области
const double Ly = 300 * units::nm;
const int Nx = 100; // количество бинов
const double carrier_density_2d = 1e16 / units::m / units::m; // m^-2 -> um^-2

struct AcousticScattering : public Scattering {
  double Da, s;  // параметры материала
  
  AcousticScattering(const Material &m, double Da_, double s_)
      : Scattering(m, 0), Da(Da_), s(s_) {}
    
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p);  // энергия частицы (с учетом анизотропии и щели)
    double T = m.get_temperature(x);
    double kT = consts::kB * T;
    // W_ac = (2 D_a² kT E)/(π ℏ³ ρ s²) * √((2 m_x m_y)/Δ(Δ+E)) * K((Δ-E)/(Δ+E))
    
    // Константа без E (E будет в rate())
    // (8π D_a² kT)/(h² ℏ ρ s²)
    double constant = 2 * Da * Da * kT / (math::pi * pow(consts::hbar, 3) * density * s * s);
    
    double Dplus = m.Delta + E;
    double Dminus = m.Delta - E;
    
    // √((2 m_x m_y)/Δ(Δ+E))
    double sqrt_term = std::sqrt((2 * m.mx * m.my) / (m.Delta * Dplus));
    
    // Параметр для эллиптического интеграла k = (Δ-E)/(Δ+E)
    double k_param = Dminus / Dplus;
    
    // Полный эллиптический интеграл первого рода K(k)
    double K = math::comp_ellint_1(k_param);
    
    return constant * E * sqrt_term * K;
  }
};

struct OpticalScattering : public Scattering {
  double constant;
  double Do;
  double phonon_energy;
  
  OpticalScattering(const Material &m, double Do_, double phonon_energy_, double signed_energy_transfer)
      : Scattering(m, signed_energy_transfer), Do(Do_), phonon_energy(phonon_energy_) {
    double omega0 = phonon_energy_ / consts::hbar;
    
    // W_op = (D_o² E) / (π ℏ² ρ ω₀) * √((2 m_x m_y)/Δ(Δ+E)) * K((Δ-E)/(Δ+E))
    // Константа без E (E будет в rate())
    constant = Do * Do / (math::pi * std::pow(consts::hbar, 2) * density * omega0);
  }

  double final_state_factor(double final_E) const {
    double Dplus = m.Delta + final_E;
    double Dminus = m.Delta - final_E;
    
    double sqrt_term = std::sqrt((2 * m.mx * m.my) / (m.Delta * Dplus));
    double k_param = Dminus / Dplus;
    double K = math::comp_ellint_1(k_param);
    
    return final_E * sqrt_term * K;
  }

  double bose(double T) const {
    return 1.0 / std::expm1(phonon_energy / (consts::kB * T));
  }
};

// Испускание оптического фонона
struct OpticalEmissionScattering : public OpticalScattering {
  OpticalEmissionScattering(const Material &m, double Do_, double phonon_energy)
      : OpticalScattering(m, Do_, phonon_energy, phonon_energy) {}
  
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p) - energy;
    if (E < m.Delta) return 0;
    return constant * final_state_factor(E) * (bose(m.get_temperature(x)) + 1);
  }
};

// Поглощение оптического фонона
struct OpticalAbsorptionScattering : public OpticalScattering {
  OpticalAbsorptionScattering(const Material &m, double Do_, double phonon_energy)
      : OpticalScattering(m, Do_, phonon_energy, -phonon_energy) {}
  
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p) - energy;
    if (E < m.Delta) return 0;
    return constant * final_state_factor(E) * bose(m.get_temperature(x));
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

void write_temperature_profile(const std::string &path,
                               const Material &material,
                               const std::vector<double> &bin_excess_energy,
                               const std::vector<uint64_t> &bin_samples,
                               double T_left,
                               double T_right) {
  std::ofstream temp_file(path);
  double T_min = 1 * units::K;
  double T_max = 3 * std::max(T_left, T_right);
  temp_file << "# x_m prescribed_T_K mean_excess_eV measured_T_K samples\n";
  for (int bin = 0; bin < material.Nx; ++bin) {
    double x_center = (bin + 0.5) * material.cell_size;
    double prescribed_T = material.get_temperature(x_center);
    double mean_excess = bin_samples[bin] > 0 ? bin_excess_energy[bin] / bin_samples[bin] : 0;
    double measured_T = bin_samples[bin] > 0
        ? material.temperature_from_mean_excess_energy(mean_excess, T_min, T_max)
        : 0;
    temp_file << x_center / units::m << " "
              << prescribed_T / units::K << " "
              << mean_excess / units::eV << " "
              << measured_T / units::K << " "
              << bin_samples[bin] << "\n";
  }
}

void write_boundary_injection_temperature(const std::string &path,
                                          const Material &material,
                                          double left_excess_energy,
                                          uint64_t left_samples,
                                          double right_excess_energy,
                                          uint64_t right_samples,
                                          double T_left,
                                          double T_right) {
  std::ofstream out(path);
  double T_min = 1 * units::K;
  double T_max = 3 * std::max(T_left, T_right);
  out << "# side prescribed_T_K mean_excess_eV measured_T_K samples\n";

  double left_mean = left_samples > 0 ? left_excess_energy / left_samples : 0;
  double left_T = left_samples > 0
      ? material.temperature_from_mean_flux_excess_energy(left_mean, T_min, T_max)
      : 0;
  out << "left "
      << T_left / units::K << " "
      << left_mean / units::eV << " "
      << left_T / units::K << " "
      << left_samples << "\n";

  double right_mean = right_samples > 0 ? right_excess_energy / right_samples : 0;
  double right_T = right_samples > 0
      ? material.temperature_from_mean_flux_excess_energy(right_mean, T_min, T_max)
      : 0;
  out << "right "
      << T_right / units::K << " "
      << right_mean / units::eV << " "
      << right_T / units::K << " "
      << right_samples << "\n";
}

int main(int argc, char const *argv[]) {
  if (argc != 8) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0] << " <ensemble size> <T_left> <T_right> <Ex> <Ey> <Bz> <all_time>\n";
    return 1;
  }

  size_t ensemble_size = parse<int>(argv[1]);
  double T_left = parse<double>(argv[2]) * units::K;
  double T_right = parse<double>(argv[3]) * units::K;
  Vec2 electric_field{parse<double>(argv[4]) * units::V / units::m, 
                      parse<double>(argv[5]) * units::V / units::m};
  double magnetic_field_z = parse<double>(argv[6]) * units::T;
  double all_time = parse<double>(argv[7]) * units::s;

  Material phosphorene{
    1.285 * consts::me, // mx (ZZ)
    0.125 * consts::me, // my (AC)
    Delta,
    Lx,
    Ly,
    Nx, 
    T_left,
    T_right
  };
  
  // Вектор механизмов рассеяния
  std::vector<Scattering *> mechanisms{
    new AcousticScattering(
      phosphorene,
      acoustic_deformation_potential,
      sound_velocity
    ),
    new OpticalEmissionScattering(
      phosphorene,
      5.67e10 * units::eV / units::m, // TO2
      56e-3 * units::eV
    ),
    new OpticalAbsorptionScattering(
      phosphorene,
      5.67e10 * units::eV / units::m, // TO2
      56e-3 * units::eV
    )
  };

  double time_step = 1e-16 * units::s;

  // print info on stdout
  std::cout << "Ensemble size:    " << ensemble_size << "\n";
  std::cout << "X length:         " << Lx / units::m << " m\n";
  std::cout << "Y length:         " << Ly / units::m << " m\n";
  std::cout << "Time step:        " << time_step / units::s << " s\n";
  std::cout << "Simulation time:  " << all_time / units::s << " s\n";
  std::cout << "Temperature left: " << T_left / units::K << " K\n";
  std::cout << "Temperature right:" << T_right / units::K << " K\n";
  std::cout << "Electric field:   " << electric_field / units::V * units::m << " V/m\n";
  std::cout << "Magnetic field:   " << magnetic_field_z / units::T << " T\n";
  std::cout << "Scattering mechanisms:\n";
  for (size_t i = 0; i < mechanisms.size(); ++i) {
    std::cout << "  " << i + 1 << ": " << *mechanisms[i] << '\n';
  }
  auto results = simulate(phosphorene,
                          mechanisms,
                          electric_field,
                          magnetic_field_z,
                          time_step,
                          all_time,
                          ensemble_size,
                          DumpFlags(DumpFlags::heat_flux | DumpFlags::particle_flux));
  // Обработка результатов
  Vec2 average_velocity;
  Vec2 average_velocity2;
  std::vector<double> scattering_rates(mechanisms.size(), 0);
  std::vector<double> total_counts(mechanisms.size(), 0);
  double avg_heat_flux = 0;
  double avg_particle_flux = 0;
  size_t flux_samples = 0;
  std::vector<double> bin_excess_energy(Nx, 0);
  std::vector<uint64_t> bin_samples(Nx, 0);
  std::vector<double> initial_bin_excess_energy(Nx, 0);
  std::vector<uint64_t> initial_bin_samples(Nx, 0);
  double injected_left_excess_energy = 0;
  double injected_right_excess_energy = 0;
  uint64_t injected_left_samples = 0;
  uint64_t injected_right_samples = 0;
  
  for (std::size_t i = 0; i < results.size(); ++i) {
    average_velocity += (results[i].average_velocity - average_velocity) / (i + 1);
    average_velocity2 +=
        (results[i].average_velocity * results[i].average_velocity - average_velocity2) / (i + 1);
        
    for (std::size_t j = 0; j < mechanisms.size(); ++j) {
      scattering_rates[j] +=
          (results[i].scattering_count[j] / all_time * units::s - scattering_rates[j]) / (i + 1);
	  total_counts[j] += results[i].scattering_count[j];
    }
    avg_heat_flux += results[i].heat_flux_sum;
    avg_particle_flux += results[i].particle_flux_sum;
    flux_samples += results[i].flux_samples;
    for (std::size_t j = 0; j < results[i].bin_excess_energy.size(); ++j) {
      bin_excess_energy[j] += results[i].bin_excess_energy[j];
      bin_samples[j] += results[i].bin_samples[j];
      initial_bin_excess_energy[j] += results[i].initial_bin_excess_energy[j];
      initial_bin_samples[j] += results[i].initial_bin_samples[j];
    }
    injected_left_excess_energy += results[i].injected_left_excess_energy;
    injected_right_excess_energy += results[i].injected_right_excess_energy;
    injected_left_samples += results[i].injected_left_samples;
    injected_right_samples += results[i].injected_right_samples;
  }
  if (flux_samples > 0) {
    avg_heat_flux /= flux_samples;
    avg_particle_flux /= flux_samples;
  }
  Vec2 std_velocity = (average_velocity2 - average_velocity * average_velocity).sqrt();
  
  // сохранение оценок теплового потока и kappa в файл
  std::filesystem::create_directories("../output");
  std::ofstream flux_file("../output/heat_flux_kappa_avg.txt");
  double gradT = ((T_right - T_left) / units::K) / (Lx / units::m);
  double heat_flux_unit_2d = units::J / units::s / units::m;
  size_t flux_windows = results.empty() ? 0 : results[0].flux_window_samples.size();
  size_t flux_stride = results.empty() ? 1 : results[0].flux_sample_stride;
  for (size_t j = 0; j < flux_windows; ++j) {
    double sum_heat_flux = 0;
    double sum_particle_flux = 0;
    uint64_t window_samples = 0;
    for (size_t i = 0; i < results.size(); ++i) {
      sum_heat_flux += results[i].heat_flux_windows[j];
      sum_particle_flux += results[i].particle_flux_windows[j];
      window_samples += results[i].flux_window_samples[j];
    }
    if (window_samples == 0) {
      continue;
    }
    double heat_flux_2d = (sum_heat_flux / window_samples) * carrier_density_2d / heat_flux_unit_2d;
    double particle_flux_2d = (sum_particle_flux / window_samples) * carrier_density_2d;
    double kappa_2d = -heat_flux_2d / gradT;
    flux_file << j * flux_stride << " " << heat_flux_2d << " " << particle_flux_2d << " " << kappa_2d << "\n";
  }
  flux_file.close();  

  write_temperature_profile("../output/initial_temperature_profile.txt",
                            phosphorene,
                            initial_bin_excess_energy,
                            initial_bin_samples,
                            T_left,
                            T_right);
  write_temperature_profile("../output/local_temperature_profile.txt",
                            phosphorene,
                            bin_excess_energy,
                            bin_samples,
                            T_left,
                            T_right);
  write_boundary_injection_temperature("../output/boundary_injection_temperature.txt",
                                       phosphorene,
                                       injected_left_excess_energy,
                                       injected_left_samples,
                                       injected_right_excess_energy,
                                       injected_right_samples,
                                       T_left,
                                       T_right);
  
  // сохранение начальных энергий частиц
//   std::ofstream energy_file("../output/initial_energy.txt");
//   for (size_t i = 0; i < results.size(); ++i) {
//     if (!results[i].energies.empty()) {
//       energy_file << results[i].energies[0] / units::eV << "\n";
//     }
//   }
//   energy_file.close();
  
  std::cout << "\n===== Results =====\n";
  std::cout << "      Directions: {ZZ, AC}\n";
  std::cout << "Average velocity: " << average_velocity << " μm/ps\n";
  std::cout << "             std: " << std_velocity << " μm/ps\n";
  double heat_flux_2d = avg_heat_flux * carrier_density_2d / heat_flux_unit_2d;
  double kappa_2d = -heat_flux_2d / gradT;
  std::cout << "Mean heat flux:   " << heat_flux_2d << " W/m (2D sheet)\n";
  std::cout << "Particle flux:    " << avg_particle_flux * carrier_density_2d << " 1/(μm·ps)\n";
  std::cout << "Kappa 2D:         " << kappa_2d << " W/K\n";
  std::cout << "Scattering rates: " << scattering_rates << " 1/s\n";
  std::cout << "           Total: " << sum(scattering_rates) << " 1/s\n";
  std::cout << "Scattering count: " << total_counts << " times\n";
  
  return 0;
}
