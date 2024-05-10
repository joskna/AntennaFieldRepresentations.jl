<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.svg" height="190">
  <img alt="" src="" height="190">
</picture>

[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joskna.github.io/AntennaFieldRepresentations.jl/dev/)

# AntennaFieldTransformations

## Introduction
This package contains various representations of the electromagnetic fields radiated by antennas.
Most importantly, the package implements the linear operator between any (parametrized) antenna representation and a predefined field sampling. 
The declared goal is to provide (efficient) methods to transfrom anyantenna representation into another to be able to choose the most suitable (most efficient) antenna representation for the problem at hand.
The purpose of this package is primarily academic and educational. At the moment and the package is under heavy development and represents by no means a finished state. Use with care!

## Featured Antenna Representations
- Collection of elctric and/or magnetic dipoles
- Vector Spherical Mode Expansion
- Far-field Pattern
- Equivalent Surface Currents (Rao-Wilton-Glisson basis)

## Featured Field Samplings
- Electric and/or magnetic field components at arbitrary (near-field) locations
- Regular spherical sampling with
  - first order probes
  - arbitrary order probes
- Irregularly sampled near-fields using probe antennas with arbitrary patterns

## Installation
The package is not registered in the Julia General Registry, yet, but a future version of this package might.
For the time being, you can install the package by switching into the package manager (type `]` in the REPL) and run
```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl
```
