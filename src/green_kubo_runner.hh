#pragma once

#include <string>
#include <vector>

#include <probable.hh>

namespace probable {

struct GreenKuboConfig {
  size_t ensemble_size;
  Vec2 electric_field;
  double magnetic_field_z;
  double all_time;
  double time_step;
  double carrier_density_2d;
  double temperature;
  std::string output_dir = "../output";
  size_t sample_stride = 100;
  size_t max_correlation_lag = 1000;
  size_t blocks = 8;
  bool pauli_time_histogram = false;
  size_t pauli_bins = 24;
  double pauli_preparation_time = 16e-9 * units::s;
  size_t pauli_preparation_passes = 12;
  size_t pauli_histories_per_pass = 64;
  std::string pauli_occupation_file;
};

struct GreenKuboRunResult {
  double kappa_2d = 0;
  double kappa_j0_2d = 0;
  double kappa_std_2d = 0;
  double kappa_j0_std_2d = 0;
  double particle_flux_2d = 0;
  double heat_flux_2d = 0;
  size_t samples = 0;
  size_t correlation_lags = 0;
  std::vector<double> scattering_rates;
  std::vector<double> scattering_counts;
  double pauli_histogram_outside_fraction = 0;
  double pauli_histogram_clipped_density_fraction = 0;
  size_t pauli_histogram_saturated_cells = 0;
};

GreenKuboRunResult run_green_kubo_simulation(const Material &material,
                                             const std::vector<Scattering *> &mechanisms,
                                             const GreenKuboConfig &config);

} // namespace probable
