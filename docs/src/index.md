# AntennaFieldRepresentations.jl

This package contains various representations of the electromagnetic fields radiated by antennas.
The declared goal is to implement the linear operator between any (parametrized) antenna representation and a predefined field sampling, and possibly its inverse to solve electromagnetic source reconstruction problems. 
To reach this goal, (efficient) methods to transform any antenna representation into another are provided. 
The possibility to switch between antenna field representations enables the user to choose the most suitable (most efficient) representation for the problem at hand.
The purpose of this package is primarily academic and educational. 
The package is under heavy development at the moment and represents by no means a finished state. Use with care!

!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` and SI units are implied everywhere.
---

## Installation

The package is not registered in the Julia General Registry, yet, but a future version of this package might.
For the time being, you can install the package by switching into the package manager (type `]` in the REPL) and run
```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl 
```

## Featured Antenna Representations
The following types of antenna representations are provided by `AntennaFieldRepresentations`:


## Featured Field Samplings
The following types of field samplings can be used with any antenna representation:
- Electric and/or magnetic field components at arbitrary (near-field) locations
- Regular spherical sampling with
  - first order probes
  - arbitrary order probes
- Irregularly sampled near-fields using probe antennas with arbitrary patterns

## Featured Operations on Antenna Representations
The main purpose of `AntennaFieldRepresentations` is to calculate the interaction between a given antenna under test (AUT) and a predefined field sampling, i.e., to calculate the received signal of a probe antenna at predefined positions (the term *position* refers to the probe's location and orientation) when the AUT is radiating. However, the functionality of the package `AntennaFieldRepresentations` is not limited to this. Furthermore, it can
- calculate the near- and far-fields of the AUT at arbitrary positions
- convert an antenna field representation into a representation of different type
- convert an antenna field representation into a representation of the same type in a different coordinate system (the new coordinate system may be rotated or translated from the original coordinate system)
- evaluate the transposed, Hermitian, and inverses of most featured linear operators
- calculate the expansion coefficients of an antenna representation given a field sampling and the received signals (i.e, solve a source reconstruction problem)
- accelerate the evaluation of the interaction operator beween suitable antenna representations and field samplings by means of the multilevel fast multipole method
