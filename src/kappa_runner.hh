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
};

struct KappaRunResult {
  Vec2 average_velocity;
  Vec2 std_velocity;
  std::vector<double> scattering_rates;
  std::vector<double> scattering_counts;
  double heat_flux_2d = 0;
  double particle_flux_2d = 0;
  double kappa_2d = 0;
};

KappaRunResult run_kappa_simulation(const Material &material,
                                    const std::vector<Scattering *> &mechanisms,
                                    const KappaRunConfig &config);

} // namespace probable
