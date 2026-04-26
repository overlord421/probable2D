#include <kappa_runner.hh>

#include <algorithm>
#include <filesystem>
#include <fstream>

namespace probable {

namespace {

template <typename T> T sum(std::vector<T> values) {
  T total{};
  for (auto &value : values) {
    total += value;
  }
  return total;
}

void write_temperature_profile(const std::string &path,
                               const Material &material,
                               const std::vector<double> &bin_excess_energy,
                               const std::vector<uint64_t> &bin_samples) {
  std::ofstream temp_file(path);
  double T_min = 1 * units::K;
  double T_max = 3 * std::max(material.T_left, material.T_right);
  temp_file << "# x_m prescribed_T_K mean_excess_eV measured_T_K samples\n";
  for (int bin = 0; bin < material.Nx; ++bin) {
    double coord_center = (bin + 0.5) * material.cell_size;
    double prescribed_T = material.get_temperature(coord_center);
    double mean_excess = bin_samples[bin] > 0 ? bin_excess_energy[bin] / bin_samples[bin] : 0;
    double measured_T = bin_samples[bin] > 0
        ? material.temperature_from_mean_excess_energy(mean_excess, T_min, T_max)
        : 0;
    temp_file << coord_center / units::m << " "
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
                                          uint64_t right_samples) {
  std::ofstream out(path);
  double T_min = 1 * units::K;
  double T_max = 3 * std::max(material.T_left, material.T_right);
  out << "# side prescribed_T_K mean_excess_eV measured_T_K samples\n";

  double left_mean = left_samples > 0 ? left_excess_energy / left_samples : 0;
  double left_T = left_samples > 0
      ? material.temperature_from_mean_flux_excess_energy(left_mean, T_min, T_max)
      : 0;
  out << "left "
      << material.T_left / units::K << " "
      << left_mean / units::eV << " "
      << left_T / units::K << " "
      << left_samples << "\n";

  double right_mean = right_samples > 0 ? right_excess_energy / right_samples : 0;
  double right_T = right_samples > 0
      ? material.temperature_from_mean_flux_excess_energy(right_mean, T_min, T_max)
      : 0;
  out << "right "
      << material.T_right / units::K << " "
      << right_mean / units::eV << " "
      << right_T / units::K << " "
      << right_samples << "\n";
}

} // namespace

KappaRunResult run_kappa_simulation(const Material &material,
                                    const std::vector<Scattering *> &mechanisms,
                                    const KappaRunConfig &config) {
  auto results = simulate(material,
                          mechanisms,
                          config.electric_field,
                          config.magnetic_field_z,
                          config.time_step,
                          config.all_time,
                          config.ensemble_size,
                          DumpFlags(DumpFlags::heat_flux | DumpFlags::particle_flux));

  KappaRunResult result;
  Vec2 average_velocity2;
  result.scattering_rates.assign(mechanisms.size(), 0);
  result.scattering_counts.assign(mechanisms.size(), 0);

  double avg_heat_flux = 0;
  double avg_particle_flux = 0;
  size_t flux_samples = 0;
  std::vector<double> bin_excess_energy(material.Nx, 0);
  std::vector<uint64_t> bin_samples(material.Nx, 0);
  std::vector<double> initial_bin_excess_energy(material.Nx, 0);
  std::vector<uint64_t> initial_bin_samples(material.Nx, 0);
  double injected_left_excess_energy = 0;
  double injected_right_excess_energy = 0;
  uint64_t injected_left_samples = 0;
  uint64_t injected_right_samples = 0;

  for (std::size_t i = 0; i < results.size(); ++i) {
    result.average_velocity += (results[i].average_velocity - result.average_velocity) / (i + 1);
    average_velocity2 +=
        (results[i].average_velocity * results[i].average_velocity - average_velocity2) / (i + 1);

    for (std::size_t j = 0; j < mechanisms.size(); ++j) {
      result.scattering_rates[j] +=
          (results[i].scattering_count[j] / config.all_time * units::s - result.scattering_rates[j]) / (i + 1);
      result.scattering_counts[j] += results[i].scattering_count[j];
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
  result.std_velocity = (average_velocity2 - result.average_velocity * result.average_velocity).sqrt();

  std::filesystem::create_directories(config.output_dir);
  std::ofstream flux_file(config.output_dir + "/heat_flux_kappa_avg.txt");
  double gradT = ((material.T_right - material.T_left) / units::K) / (material.axis_length() / units::m);
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
    double heat_flux_2d = (sum_heat_flux / window_samples) * config.carrier_density_2d / heat_flux_unit_2d;
    double particle_flux_2d = (sum_particle_flux / window_samples) * config.carrier_density_2d;
    double kappa_2d = -heat_flux_2d / gradT;
    flux_file << j * flux_stride << " " << heat_flux_2d << " " << particle_flux_2d << " " << kappa_2d << "\n";
  }

  write_temperature_profile(config.output_dir + "/initial_temperature_profile.txt",
                            material,
                            initial_bin_excess_energy,
                            initial_bin_samples);
  write_temperature_profile(config.output_dir + "/local_temperature_profile.txt",
                            material,
                            bin_excess_energy,
                            bin_samples);
  write_boundary_injection_temperature(config.output_dir + "/boundary_injection_temperature.txt",
                                       material,
                                       injected_left_excess_energy,
                                       injected_left_samples,
                                       injected_right_excess_energy,
                                       injected_right_samples);

  result.heat_flux_2d = avg_heat_flux * config.carrier_density_2d / heat_flux_unit_2d;
  result.particle_flux_2d = avg_particle_flux * config.carrier_density_2d;
  result.kappa_2d = -result.heat_flux_2d / gradT;

  return result;
}

} // namespace probable
