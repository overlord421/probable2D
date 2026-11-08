#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <probable.hh>

namespace probable {

// Uniform physical momentum cells. Visits are weighted by elapsed time,
// including dwell intervals without an accepted scattering event.
class PauliTimeHistogram {
 public:
  PauliTimeHistogram(size_t bins, double px_limit, double py_limit,
                     double density, double spin)
      : bins_(bins), px_limit_(px_limit), py_limit_(py_limit), density_(density) {
    if (bins < 2 || bins > 512 || !(px_limit > 0) || !(py_limit > 0) ||
        !(density > 0) || !(spin > 0))
      throw std::invalid_argument("Invalid Pauli momentum histogram geometry");
    dp_x_ = 2 * px_limit / bins;
    dp_y_ = 2 * py_limit / bins;
    states_per_area_ = spin * dp_x_ * dp_y_ / std::pow(2 * math::pi * consts::hbar, 2);
    residence_.assign(bins * bins, 0);
    occupancy_.assign(bins * bins, 0);
  }

  void record(const Vec2 &p, double duration) {
    if (!(duration > 0) || !std::isfinite(duration))
      throw std::invalid_argument("Pauli visit duration must be finite and positive");
    total_time_ += duration;
    size_t cell = index(p);
    if (cell < residence_.size()) residence_[cell] += duration;
    else outside_time_ += duration;
  }

  // Homogeneous zero-field equilibrium is even in both momentum components.
  // Share each time observation between the four equivalent cells, keeping
  // its total duration unchanged. This suppresses artificial odd drift.
  void record_reflected(const Vec2 &p, double duration) {
    if (!(duration > 0) || !std::isfinite(duration))
      throw std::invalid_argument("Pauli visit duration must be finite and positive");
    total_time_ += duration;
    size_t cell = index(p);
    if (cell >= residence_.size()) { outside_time_ += duration; return; }
    size_t x = cell / bins_, y = cell % bins_;
    for (size_t ix : {x, bins_ - 1 - x})
      for (size_t iy : {y, bins_ - 1 - y})
        residence_[ix * bins_ + iy] += duration / 4;
  }

  // Updating occurs only between preparation passes. Production reads a
  // frozen map; observations cannot silently change its acceptance rule.
  void initialize(size_t cell, double value) {
    if (!(value >= 0 && value < 1)) throw std::invalid_argument("Invalid initial occupancy");
    occupancy_.at(cell) = value;
  }
  double update(double relaxation = 1, bool vacancy_damping = false) {
    if (!(total_time_ > 0)) throw std::invalid_argument("Empty Pauli time histogram");
    if (!(relaxation > 0 && relaxation <= 1)) throw std::invalid_argument("Invalid map relaxation");
    double change = 0;
    for (size_t i = 0; i < occupancy_.size(); ++i) {
      double value = std::clamp(raw_occupation(i), 0.0, 1.0);
      double fraction = relaxation * (vacancy_damping ? 1 - occupancy_[i] : 1);
      double next = occupancy_[i] + fraction * (value - occupancy_[i]);
      change = std::max(change, std::fabs(next - occupancy_[i]));
      occupancy_[i] = next;
    }
    return change;
  }

  double occupation(const Vec2 &p) const {
    size_t cell = index(p);
    return cell < occupancy_.size() ? occupancy_[cell] : 0;
  }
  double raw_occupation(size_t i) const {
    return total_time_ > 0 ? density_ * residence_.at(i) / total_time_ / states_per_area_ : 0;
  }
  Vec2 center(size_t i) const {
    return {-px_limit_ + (i / bins_ + .5) * dp_x_,
            -py_limit_ + (i % bins_ + .5) * dp_y_};
  }
  double residence(size_t i) const { return residence_.at(i); }
  double frozen_occupation(size_t i) const { return occupancy_.at(i); }
  double total_time() const { return total_time_; }
  double outside_fraction() const { return total_time_ > 0 ? outside_time_ / total_time_ : 0; }
  double states_per_area() const { return states_per_area_; }
  double dp_x() const { return dp_x_; }
  double dp_y() const { return dp_y_; }
  double clipped_density_fraction() const {
    double mass = 0;
    for (double f : occupancy_) mass += f * states_per_area_;
    return mass / density_;
  }
  size_t cells() const { return residence_.size(); }
  size_t saturated_cells() const {
    size_t count = 0;
    for (size_t i = 0; i < cells(); ++i) if (raw_occupation(i) >= 1) ++count;
    return count;
  }

 private:
  size_t index(const Vec2 &p) const {
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || p.x < -px_limit_ ||
        p.x >= px_limit_ || p.y < -py_limit_ || p.y >= py_limit_) return residence_.size();
    size_t ix = static_cast<size_t>((p.x + px_limit_) / dp_x_);
    size_t iy = static_cast<size_t>((p.y + py_limit_) / dp_y_);
    return ix < bins_ && iy < bins_ ? ix * bins_ + iy : residence_.size();
  }
  size_t bins_;
  double px_limit_, py_limit_, density_, dp_x_, dp_y_, states_per_area_;
  double total_time_ = 0, outside_time_ = 0;
  std::vector<double> residence_, occupancy_;
};

} // namespace probable
