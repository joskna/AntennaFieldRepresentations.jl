# AntennaFieldRepresentations.jl

This package contains various representations of the electromagnetic fields radiated by antennas.
`AntennaFieldRepresentations.jl` aims to implement the linear operator between any of the featured antenna representation and a predefined field sampling.
By providing the possibility to switch between various antenna field representations, `AntennaFieldRepresentations.jl` enables the user to choose the most suitable (most efficient) representation for the problem at hand.
The purpose of this package is primarily academic and educational. 
The package is under heavy development at the moment and represents by no means a finished state. Use with care!

!!! note
    The package `AntennaFieldRepresentations.jl` deals with representations of the electromagnetic fields of an antenna at a single frequency.
    
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` is implicitly assumed for the relevant variables and SI units are implied everywhere.
---

## Installation

The package is not registered in the Julia General Registry, yet, but a future version of this package might.
For the time being, you can install the package by switching into the package manager (type `]` in the REPL) and run
```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl 
```

## Featured Antenna Representations
The following types of antenna representations are provided by `AntennaFieldRepresentations`:
- Collections of short electric and magnetic dipoles
- Spherical Vector-Wave Expansions
- Plane-WaveExpansions
- Equivalent Surface-Currents


## Featured Field Samplings
The following types of field samplings can be used with any antenna representation:
- Electric and/or magnetic field components at arbitrary (near-field) locations
- Regular spherical sampling with
  - first order probes
  - arbitrary order probes
- Irregularly sampled near-fields using probe antennas with arbitrary receiving patterns

## Featured Operations on Antenna Representations
The main purpose of `AntennaFieldRepresentations.jl` is to calculate the interaction between a given antenna under test (AUT) and a predefined field sampling, i.e., to calculate the received signal of probe antennas at predefined positions when the AUT is radiating.

However, the functionality of `AntennaFieldRepresentations.jl` is not limited to this. 
With `AntennaFieldRepresentations.jl` the user can (besides other useful things)
- calculate the near- and far-fields of the AUT at arbitrary positions
- convert an antenna field representation into a representation of different type
- convert an antenna field representation into a representation of the same type in a different coordinate system (the new coordinate system may be rotated or translated from the original coordinate system)
- evaluate the transposed, Hermitian, and inverses of most featured linear operators
- calculate the expansion coefficients of an antenna representation given a field sampling and the received signals (i.e, solve a source reconstruction problem)
- accelerate the evaluation of the interaction operator beween suitable antenna representations and field samplings by means of the multilevel fast multipole method

## Known Issues

- Resampling works only if for every sampling point in the target `SphereSamplingStrategy` the point on the opposite side of the sphere is also included in the sampling.
- Resampling only works if the original `SphereSamplingStrategy` is a `GaussLegendreθRegularϕSampling` with an even number of `ϕ`-samples
