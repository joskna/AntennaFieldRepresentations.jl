# AntennaFieldRepresentations.jl

This package contains various representations of the electromagnetic fields radiated by antennas.
Most importantly, the package implements the linear operator between any (parametrized) antenna representation and a predefined field sampling. 
The declared goal is to provide (efficient) methods to transfrom anyantenna representation into another to be able to choose the most suitable (most efficient) antenna representation for the problem at hand.
The purpose of this package is primarily academic and educational. At the moment and the package is under heavy development and represents by no means a finished state. Use with care!

!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` and SI units are used everywhere.
---

## Installation

Install `AntennaFieldRepresentations` by cloning the package to your ~\\.julia\dev folder and registering the package in the package manager (enter `]` at the julia REPL) and type 

```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl 
```
to register the package to your `Julia` environment. Make sure that the required related packages are available.

## Featured Antenna Representations
The following types of antenna representations are provided by `AntennaFieldRepresentations`:
- Collection of elctric and/or magnetic dipoles
- Vector Spherical Mode Expansion
- Far-field Pattern
- Equivalent Surface Currents (Rao-Wilton-Glisson basis)

## Featured Field Samplings
The following types of field samplings can be used with any antenna representation:
- Electric and/or magnetic field components at arbitrary (near-field) locations
- Regular spherical sampling with
  - first order probes
  - arbitrary order probes
- Irregularly sampled near-fields using probe antennas with arbitrary patterns

## Operations on Antenna Representations
The main purpose of `AntennaFieldRepresentations` is to calculate the interaction between a given antenna under test (AUT) and a predefined field sampling, i.e., to calculate the received signal of a probe antenna at predefined positions (the term *position* refers to the probe's location and orientation) when the AUT is radiating. However, the functionality of the package `AntennaFieldRepresentations` is not limited to this. Furthermore, it can
- calculate the near- and far-fields of the AUT at arbitrary positions
- convert an antenna field representation into a representation of different type
- convert an antenna field representation into a representation of the same type in a different coordinate system (the new coordinate system may be rotated or translated from the original coordinate system)
- calculate the expansion coefficients of an antenna representation given a field sampling and the received signals (i.e, solve a source reconstruction problem)


### Planned Features
- ⏳[^1] **Scalar spherical mode expansions** 
- ⏳[^1] **Inhomogeneous plane wave expansions**
- ⏳[^1] **Cylindrical vector modes**
- ⏳[^1] **Spheroidal wave functions**
- ⏳[^1] **All equivalent surface current types provided by** `Beast.jl`
- ⏳[^1] **Efficient implementations of Love currents and Caldéron projectors**
- ⏳[^1] Lossless dielectrics (``\varepsilon_\mathrm{r},\mu_{\mathrm{r}} \in \mathbb{R}``)
- ⏳[^1] Lossy dielectrics (``\varepsilon_\mathrm{r},\mu_{\mathrm{r}} \in \mathbb{C}``)


## [Footnotes](@id refs)
[^1]: ⏳ Planned but not yet implemented