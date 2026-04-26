#pragma once

#include <probable.hh>

namespace probable {

struct AcousticScattering : public Scattering {
  double density, Da, s;
  
  AcousticScattering(const Material &m, double density_, double Da_, double s_)
      : Scattering(m, 0), density(density_), Da(Da_), s(s_) {}
    
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p);
    double T = m.get_temperature(x);
    double kT = consts::kB * T;
    double constant = 2 * Da * Da * kT / (math::pi * pow(consts::hbar, 3) * density * s * s);
    
    double Dplus = m.Delta + E;
    double Dminus = m.Delta - E;
    double sqrt_term = std::sqrt((2 * m.mx * m.my) / (m.Delta * Dplus));
    double k_param = Dminus / Dplus;
    double K = math::comp_ellint_1(k_param);
    
    return constant * E * sqrt_term * K;
  }
};

struct OpticalScattering : public Scattering {
  double density;
  double constant;
  double Do;
  double phonon_energy;
  
  OpticalScattering(const Material &m,
                    double density_,
                    double Do_,
                    double phonon_energy_,
                    double signed_energy_transfer)
      : Scattering(m, signed_energy_transfer),
        density(density_),
        Do(Do_),
        phonon_energy(phonon_energy_) {
    double omega0 = phonon_energy_ / consts::hbar;
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

struct OpticalEmissionScattering : public OpticalScattering {
  OpticalEmissionScattering(const Material &m, double density, double Do, double phonon_energy)
      : OpticalScattering(m, density, Do, phonon_energy, phonon_energy) {}
  
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p) - energy;
    if (E < m.Delta) return 0;
    return constant * final_state_factor(E) * (bose(m.get_temperature(x)) + 1);
  }
};

struct OpticalAbsorptionScattering : public OpticalScattering {
  OpticalAbsorptionScattering(const Material &m, double density, double Do, double phonon_energy)
      : OpticalScattering(m, density, Do, phonon_energy, -phonon_energy) {}
  
  double rate(const Vec2 &p, double x) const override {
    double E = m.energy(p) - energy;
    if (E < m.Delta) return 0;
    return constant * final_state_factor(E) * bose(m.get_temperature(x));
  }
};

} // namespace probable
