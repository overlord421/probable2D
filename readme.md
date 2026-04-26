<h1 align=center>probable2D</h1>

<div align=center><em>Library for Monte Carlo simulations of electronic transport in 2D solids</em></div>

<br>

Build:  
`g++ -fopenmp -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc kappa_runner.cc -lm -lstdc++`

Portable build without OpenMP:
`g++ -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc kappa_runner.cc -lm -lstdc++`

Run with optional thermal axis and open-circuit Seebeck field fitting:
`./phosphorene <ensemble size> <T_left> <T_right> <Ex> <Ey> <Bz> <all_time> [axis] [--seebeck]`

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

Reusable code for new materials:
- `probable.hh` / `probable.cc`: EMC core, geometry, particle motion, boundary
  thermostats, and generic `Scattering` interface.
- `kappa_runner.hh` / `kappa_runner.cc`: common kappa workflow, ensemble
  averaging, flux history, temperature profiles, output files, and optional
  Seebeck field fitting for `J=0`.
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
time windows instead of storing every particle step.
Optical phonon scattering includes both emission and absorption with Bose
factors, and boundary injection diagnostics use flux-weighted temperature
reconstruction.
