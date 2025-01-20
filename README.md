<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.svg" height="190">
  <img alt="" src="" height="190">
</picture>

[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joskna.github.io/AntennaFieldRepresentations.jl/dev/)
[![Build Status](https://github.com/joskna/AntennFieldRepresentations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joskna/AntennaFieldRepresentations.jl/actions/workflows/CI.yml?query=branch%3Amain)

# AntennaFieldRepresentations

## Introduction
This package contains various representations of the electromagnetic fields radiated by antennas.
The declared goal is to implement the linear operator between any (parametrized) antenna representation and a predefined field sampling, and possibly its inverse to solve electromagnetic source reconstruction problems. 
To reach this goal, (efficient) methods to transform any antenna representation into another are provided. 
The possibility to switch between antenna field representations enables the user to choose the most suitable (most efficient) representation for the problem at hand.
The purpose of this package is primarily academic and educational. 
The package is under heavy development at the moment and represents by no means a finished state. Use with care!

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
