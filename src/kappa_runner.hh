#pragma once

#include <string>
#include <vector>

#include <probable.hh>

namespace probable {

struct KappaRunConfig {
  size_t ensemble_size;
  Vec2 electric_field;
  double magnetic_field_z;
  double all_time;
  double time_step;
  double carrier_density_2d;
  std::string output_dir = "../output";
  bool tune_seebeck_field = false;
  int seebeck_iterations = 4;
  double seebeck_initial_step = 1e4 * units::V / units::m;
  double seebeck_tolerance = 1e-2;
};

struct KappaRunResult {
  Vec2 average_velocity;
  Vec2 std_velocity;
  std::vector<double> scattering_rates;
  std::vector<double> scattering_counts;
  double heat_flux_2d = 0;
  double particle_flux_2d = 0;
  double kappa_2d = 0;
  double kappa_std_2d = 0;
  Vec2 electric_field;
  double seebeck_field_axis = 0;
};

KappaRunResult run_kappa_simulation(const Material &material,
                                    const std::vector<Scattering *> &mechanisms,
                                    const KappaRunConfig &config);

KappaRunResult run_open_circuit_kappa_simulation(const Material &material,
                                                 const std::vector<Scattering *> &mechanisms,
                                                 const KappaRunConfig &config);

} // namespace probable
