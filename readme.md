<h1 align=center>probable2D</h1>

<div align=center><em>Library for Monte Carlo simulations of electronic transport in 2D solids</em></div>

<br>

Build:  
`g++ -fopenmp -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc -lm -lstdc++`

Portable build without OpenMP:
`g++ -O2 -I. -std=c++17 -o phosphorene phosphorene.cc probable.cc -lm -lstdc++`

The reworked kappa branch also writes `output/initial_temperature_profile.txt`
and `output/local_temperature_profile.txt` with temperature profiles reconstructed
from binned carrier energy.
`output/boundary_injection_temperature.txt` reports the reconstructed
temperature of particles injected from the left and right thermostats.
Temperature and flux tallies skip the first 20% of time steps as burn-in.
Flux history in `output/heat_flux_kappa_avg.txt` is accumulated in at most 1000
time windows instead of storing every particle step.
