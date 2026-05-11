#include <iostream>
#include <sstream>
#include <string>

#include <gapped2d_scattering.hh>
#include <green_kubo_runner.hh>
#include <kappa_runner.hh>
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
  if (argc < 8) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0]
              << " <ensemble size> <T_left> <T_right> <Ex> <Ey> <Bz> <all_time>"
              << " [axis=x|y] [--seebeck|--green-kubo] [--fermi-dirac]"
              << " [--density-cm2=value] [--optical-energy-mev=value]"
              << " [--acoustic-da-ev=value] [--optical-do-ev-m=value]\n";
    return 1;
  }

  size_t ensemble_size = parse<int>(argv[1]);
  double T_left = parse<double>(argv[2]) * units::K;
  double T_right = parse<double>(argv[3]) * units::K;
  Vec2 electric_field{parse<double>(argv[4]) * units::V / units::m, 
                      parse<double>(argv[5]) * units::V / units::m};
  double magnetic_field_z = parse<double>(argv[6]) * units::T;
  double all_time = parse<double>(argv[7]) * units::s;
  char thermal_axis = 'x';
  bool tune_seebeck_field = false;
  bool green_kubo = false;
  bool use_fermi_dirac = false;
  double run_carrier_density_2d = carrier_density_2d;
  double run_acoustic_deformation_potential = acoustic_deformation_potential;
  double run_optical_deformation_potential = 5.67e10 * units::eV / units::m;
  double run_optical_phonon_energy = 56e-3 * units::eV;
  for (int i = 8; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "y" || arg == "Y" || arg == "axis=y") {
      thermal_axis = 'y';
    } else if (arg == "x" || arg == "X" || arg == "axis=x") {
      thermal_axis = 'x';
    } else if (arg == "--seebeck" || arg == "seebeck") {
      tune_seebeck_field = true;
    } else if (arg == "--green-kubo" || arg == "--gk" || arg == "green-kubo" || arg == "gk") {
      green_kubo = true;
    } else if (arg == "--fermi-dirac" || arg == "--fd" || arg == "fermi-dirac" || arg == "fd") {
      use_fermi_dirac = true;
    } else if (arg.rfind("--density-cm2=", 0) == 0 || arg.rfind("--carrier-density-cm2=", 0) == 0) {
      std::string value = arg.substr(arg.find('=') + 1);
      run_carrier_density_2d = parse<double>(value) * 1e4 / units::m / units::m;
    } else if (arg.rfind("--optical-energy-mev=", 0) == 0) {
      std::string value = arg.substr(arg.find('=') + 1);
      run_optical_phonon_energy = parse<double>(value) * 1e-3 * units::eV;
    } else if (arg.rfind("--acoustic-da-ev=", 0) == 0) {
      std::string value = arg.substr(arg.find('=') + 1);
      run_acoustic_deformation_potential = parse<double>(value) * units::eV;
    } else if (arg.rfind("--optical-do-ev-m=", 0) == 0) {
      std::string value = arg.substr(arg.find('=') + 1);
      run_optical_deformation_potential = parse<double>(value) * units::eV / units::m;
    } else {
      std::cout << "Unknown optional argument: " << arg << "\n";
      return 1;
    }
  }
  if (green_kubo && tune_seebeck_field) {
    std::cout << "--green-kubo and --seebeck are separate kappa estimators; choose one mode\n";
    return 1;
  }

  double equilibrium_temperature = 0.5 * (T_left + T_right);
  double material_T_left = green_kubo ? equilibrium_temperature : T_left;
  double material_T_right = green_kubo ? equilibrium_temperature : T_right;

  Material phosphorene{
    1.285 * consts::me, // mx (ZZ)
    0.125 * consts::me, // my (AC)
    Delta,
    Lx,
    Ly,
    Nx, 
    material_T_left,
    material_T_right,
    thermal_axis,
    use_fermi_dirac,
    run_carrier_density_2d
  };
  
  // Вектор механизмов рассеяния
  std::vector<Scattering *> mechanisms{
    new AcousticScattering(
      phosphorene,
      density,
      run_acoustic_deformation_potential,
      sound_velocity
    ),
    new OpticalEmissionScattering(
      phosphorene,
      density,
      run_optical_deformation_potential, // TO2
      run_optical_phonon_energy
    ),
    new OpticalAbsorptionScattering(
      phosphorene,
      density,
      run_optical_deformation_potential, // TO2
      run_optical_phonon_energy
    )
  };

  double time_step = 1e-16 * units::s;

  // print info on stdout
  std::cout << "Ensemble size:    " << ensemble_size << "\n";
  std::cout << "X length:         " << Lx / units::m << " m\n";
  std::cout << "Y length:         " << Ly / units::m << " m\n";
  std::cout << "Thermal axis:     " << phosphorene.thermal_axis << "\n";
  std::cout << "Time step:        " << time_step / units::s << " s\n";
  std::cout << "Simulation time:  " << all_time / units::s << " s\n";
  std::cout << "Temperature left: " << T_left / units::K << " K\n";
  std::cout << "Temperature right:" << T_right / units::K << " K\n";
  if (green_kubo) {
    std::cout << "Equilibrium T:    " << equilibrium_temperature / units::K << " K\n";
  }
  std::cout << "Carrier density:  " << run_carrier_density_2d * units::m * units::m / 1e4 << " cm^-2\n";
  std::cout << "Acoustic Da:      " << run_acoustic_deformation_potential / units::eV << " eV\n";
  std::cout << "Optical Do:       " << run_optical_deformation_potential / units::eV * units::m << " eV/m\n";
  std::cout << "Optical phonon:   " << run_optical_phonon_energy / units::eV * 1e3 << " meV\n";
  std::cout << "Electric field:   " << electric_field / units::V * units::m << " V/m\n";
  std::cout << "Seebeck fitting:  " << (tune_seebeck_field ? "on" : "off") << "\n";
  std::cout << "Green-Kubo mode:  " << (green_kubo ? "on" : "off") << "\n";
  std::cout << "Statistics:       " << (use_fermi_dirac ? "Fermi-Dirac" : "Boltzmann") << "\n";
  if (use_fermi_dirac) {
    std::cout << "Chemical mu:      " << phosphorene.chemical_potential(green_kubo ? equilibrium_temperature : 0.5 * (T_left + T_right)) / units::eV << " eV\n";
  }
  std::cout << "Magnetic field:   " << magnetic_field_z / units::T << " T\n";
  std::cout << "Scattering mechanisms:\n";
  for (size_t i = 0; i < mechanisms.size(); ++i) {
    std::cout << "  " << i + 1 << ": " << *mechanisms[i] << '\n';
  }
  if (green_kubo) {
    GreenKuboConfig config{
      ensemble_size,
      electric_field,
      magnetic_field_z,
      all_time,
      time_step,
      run_carrier_density_2d,
      equilibrium_temperature,
      "../output"
    };
    GreenKuboRunResult result = run_green_kubo_simulation(phosphorene, mechanisms, config);

    std::cout << "\n===== Green-Kubo Results =====\n";
    std::cout << "      Directions: {ZZ, AC}\n";
    std::cout << "Correlation samples: " << result.samples << "\n";
    std::cout << "Correlation lags:    " << result.correlation_lags << "\n";
    std::cout << "Mean heat flux:      " << result.heat_flux_2d << " W/m (2D sheet)\n";
    std::cout << "Particle flux:       " << result.particle_flux_2d << " 1/(μm·ps)\n";
    std::cout << "GK Kappa 2D:         " << result.kappa_2d << " W/K\n";
    std::cout << "GK Kappa 2D std:     " << result.kappa_std_2d << " W/K\n";
    std::cout << "GK Kappa J=0 2D:     " << result.kappa_j0_2d << " W/K\n";
    std::cout << "GK Kappa J=0 std:    " << result.kappa_j0_std_2d << " W/K\n";
    std::cout << "Scattering rates:    " << result.scattering_rates << " 1/s\n";
    std::cout << "           Total:    " << sum(result.scattering_rates) << " 1/s\n";
    std::cout << "Scattering count:    " << result.scattering_counts << " times\n";
    return 0;
  }

  KappaRunConfig config{
    ensemble_size,
    electric_field,
    magnetic_field_z,
    all_time,
    time_step,
    run_carrier_density_2d,
    "../output"
  };
  config.tune_seebeck_field = tune_seebeck_field;
  KappaRunResult result = run_open_circuit_kappa_simulation(phosphorene, mechanisms, config);

  std::cout << "\n===== Results =====\n";
  std::cout << "      Directions: {ZZ, AC}\n";
  std::cout << "Electric field:   " << result.electric_field / units::V * units::m << " V/m\n";
  std::cout << "Seebeck axis E:   " << result.seebeck_field_axis / units::V * units::m << " V/m\n";
  std::cout << "Average velocity: " << result.average_velocity << " μm/ps\n";
  std::cout << "             std: " << result.std_velocity << " μm/ps\n";
  std::cout << "Mean heat flux:   " << result.heat_flux_2d << " W/m (2D sheet)\n";
  std::cout << "Particle flux:    " << result.particle_flux_2d << " 1/(μm·ps)\n";
  std::cout << "Kappa 2D:         " << result.kappa_2d << " W/K\n";
  std::cout << "Kappa 2D std:     " << result.kappa_std_2d << " W/K\n";
  std::cout << "Scattering rates: " << result.scattering_rates << " 1/s\n";
  std::cout << "           Total: " << sum(result.scattering_rates) << " 1/s\n";
  std::cout << "Scattering count: " << result.scattering_counts << " times\n";
  
  return 0;
}
