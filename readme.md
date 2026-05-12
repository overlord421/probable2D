<h1 align=center>probable2D</h1>

<div align=center><em>Library for Monte Carlo simulations of electronic transport in 2D solids</em></div>

<br>

Build:  
`g++ -fopenmp -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc kappa_runner.cc green_kubo_runner.cc -lm -lstdc++`

Portable build without OpenMP:
`g++ -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc kappa_runner.cc green_kubo_runner.cc -lm -lstdc++`

Run with optional thermal axis, open-circuit Seebeck field fitting, or
equilibrium Green-Kubo kappa:
`./phosphorene <ensemble size> <T_left> <T_right> <Ex> <Ey> <Bz> <all_time> [axis] [--seebeck|--green-kubo] [--fermi-dirac]`

Optional material overrides for parameter sweeps:
`--density-cm2=... --optical-energy-mev=... --acoustic-da-ev=... --optical-do-ev-m=... --sound-velocity-ms=...`

For `axis=x`, thermostats are placed at `x=0` and `x=Lx`, the temperature
profile and heat flux are measured along `x`, and `y` is periodic. For
`axis=y`, thermostats are placed at `y=0` and `y=Ly`, the profile and heat flux
are measured along `y`, and `x` is periodic.

With `--seebeck`, the code reruns the simulation for several trial electric
fields along the selected thermal axis and uses a secant update to make the
mean particle flux close to zero. The reported kappa then corresponds to the
open-circuit condition `J=0` for the chosen finite sample. Trial results are
saved in `output/seebeck_field_fit.txt`; the usual temperature and flux output
files are written for the final selected field.

With `--green-kubo`, the code runs an equilibrium simulation at
`T0 = (T_left + T_right) / 2` with periodic boundaries in both directions.
It writes `output/gk_heat_current_correlation_x.txt` or `_y.txt` with the
`QQ`, `QJ`, and `JJ` current correlation functions,
`output/gk_kappa_running_x.txt` or `_y.txt` with the running Green-Kubo
integrals, and `output/gk_kappa_blocks_x.txt` or `_y.txt` with block estimates.
The reported `GK Kappa J=0 2D` uses the correlation correction
`I_QQ - I_QJ^2 / I_JJ`.

With `--fermi-dirac`, initial particles and boundary-injected particles are
sampled from Fermi-Dirac statistics at the configured sheet density. Scattering
events include Pauli blocking through the final-state factor `1 - f(E_final)`.
The chemical potential is found numerically from the actual 2D dispersion
`E(p)` and the target sheet density.
Without this flag, the original Boltzmann sampling and unblocked scattering are
used.

Reusable code for new materials:
- `probable.hh` / `probable.cc`: EMC core, geometry, particle motion, boundary
  thermostats, and generic `Scattering` interface.
- `kappa_runner.hh` / `kappa_runner.cc`: common kappa workflow, ensemble
  averaging, flux history, temperature profiles, output files, and optional
  Seebeck field fitting for `J=0`.
- `green_kubo_runner.hh` / `green_kubo_runner.cc`: equilibrium Green-Kubo
  kappa workflow based on heat-current and particle-current correlations.
- `gapped2d_scattering.hh`: reusable acoustic and optical phonon mechanisms for
  the current anisotropic gapped 2D dispersion. A new material can reuse these
  classes with different constants, or provide its own `Scattering` subclasses.
- `phosphorene.cc`: phosphorene-specific constants, material construction, and
  mechanism list.

The reworked kappa branch also writes `output/initial_temperature_profile.txt`
and `output/local_temperature_profile.txt` with temperature profiles reconstructed
from binned carrier energy.
`output/boundary_injection_temperature.txt` reports the reconstructed
temperature of particles injected from the left and right thermostats.
Temperature and flux tallies skip the first 20% of time steps as burn-in.
Flux history in `output/heat_flux_kappa_avg.txt` is accumulated in at most 1000
time windows instead of storing every particle step. Its columns are step,
mean 2D heat flux, mean particle flux, mean 2D kappa, and the standard
deviation of 2D kappa across ensemble trajectories.
Optical phonon scattering includes both emission and absorption with Bose
factors, and boundary injection diagnostics use flux-weighted temperature
reconstruction.
