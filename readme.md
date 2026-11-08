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

The default scattering set is acoustic deformation-potential scattering plus
two intraband Gamma optical modes, each with emission and absorption:

| Mode | Energy (meV) | Absolute deformation potential (eV/m) |
| --- | ---: | ---: |
| 7 | 41.786546532 | 5.802009891812528e10 |
| 11 | 54.464547517 | 5.598752616903092e10 |

The optical constants are fixed-response band-11 elements from QE job 8576;
energies use the fixed-occupation crystal-ASR frequencies from job 8568.
Both constants agree with converged frozen-phonon slopes within 0.11%.
This is the selected two-mode approximation to intraband scattering;
finite-q and interband transport are outside this parameterization.
The acoustic default is 4.9 eV and remains independently configurable.
All nine optical-mode results are retained in
`qe_phosphorene_gamma/dfpt_fixed_dvscf_tetraopt_run_8575/`.

The same two-mode set is saved in `materials/phosphorene_gamma_optical_modes.txt`.
A mode file replaces the default set. When running from `src/`, use:
`--optical-modes-file=../materials/phosphorene_gamma_optical_modes.txt`

Explicit `--optical-energy-mev` or `--optical-do-ev-m` selects the legacy
single-mode override; its unspecified value retains 56 meV or 5.67e10 eV/m,
respectively. An explicit mode file takes precedence over these overrides.

The optical modes file format is:
`name energy_meV deformation_potential_eV_m`

For each valid mode the code creates both optical emission and absorption
scattering channels. Modes with non-positive energy or deformation potential
are skipped.

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

Default working algorithm: the estimator of article run 8513, with
exact conditional energy-surface sampling and the previous tagged-electron
Green-Kubo estimator. On 2026-10-06 the user requested postponing the new
FD response algorithm; its source, checks and figures remain archived in
`transport_runs/article_linear_fd_N8_20261006_863743c/`. Those results are
not the selected working baseline. The optical defaults now use the two-mode
QE parameterization above; theoretical estimator validation remains open.

For homogeneous zero-field FD Green-Kubo, `--pauli=time-histogram` enables
occupation estimated from preparation trajectories sampled at fixed time
intervals. The histories run sequentially, independently of production
ensemble size. FD initializes particles and the first preparation pass;
later passes use a frozen empirical momentum map, with vacancy-damped
updates between passes. Production uses the final frozen map. The old
correlation formulas, carrier-density normalization and scattering rates
are unchanged. The default remains `--pauli=equilibrium`.

Preparation controls: `--pauli-bins=24`, `--pauli-prep-ns=16` (total time),
`--pauli-passes=12`, `--pauli-histories=64` (sequential histories per pass).
The output includes `pauli_preparation_x.tsv` / `_y.tsv` and
`pauli_occupation_x.tsv` / `_y.tsv`, including raw and used occupations,
FD cell-average reference, map changes, overflow and density diagnostics.
Unresolved preparation stops before measuring correlations. Grid/time
convergence remains necessary for quantitative use; these preparation
checks alone do not establish transport accuracy. See
`transport_runs/pauli_time_histogram_20261006/METHOD.md`.

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

For zero-field preparation, observations are shared between the four momentum
reflection partners without changing their total dwell time. Bare rates are
cached between accepted collisions, retaining the original step hazards.
Both production axes can reuse the same checked map via
`--pauli-map=/absolute/path/pauli_occupation_x.tsv --pauli-bins=N`. Loading
checks cell geometry, FD reference/material, density/state capacity, overflow
and reflection symmetry; it skips preparation and freezes the supplied map.

## Article driver for linearized Fermi–Dirac response

Build and run from the repository root.

```sh
run_dir=transport_runs/mc_article_density_20261020
c++ -O2 -std=c++17 -Isrc "$run_dir/optical_response.cc" src/probable.cc \
  -o "$run_dir/optical_response"
OPTICAL_RESPONSE_SEED=5102400 "$run_dir/optical_response" 2e13 131072 "$run_dir/modes_7_11.txt"
```

The arguments are sheet density in cm⁻², number of histories, and the mode
file. The JSON output includes the correlation integrals from both complete
event trajectories and conditional averaging. The article uses the
`conditional_moments` integrals. Setting `OPTICAL_RESPONSE_SEED` reproduces the
seed used by the archived driver. The generator uses a random seed by default.

Run `python3 transport_runs/mc_article_density_20261020/run_density.py` from a
fresh checkout to generate the five density series with their recorded seeds.
The script writes results under the driver's directory and requires new
case directories.
