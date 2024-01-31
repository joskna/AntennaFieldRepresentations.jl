<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.svg" height="190">
  <img alt="" src="" height="190">
</picture>

# AntennaFieldTransformations

## Introduction
This package contains various representations of electromagnetic fields and their interactions with antennas. 
The goal provide (efficient) methods to transfrom any representation into another which might be more suitable for the problem at hand.
The purpose of this package is primarily academic and educational at the moment and the package is primarily under heavy at the moment and represents by no means a finished state. Use with care!

## Installation
At the moment, the package is not registered in the Julia General Registry, but a future version of this package might.
For the time being, you can install the package by switching into the package manager (type `]` in the REPL) and run
```
add https://github.com/joskna/AntennaFieldRepresentations.jl
```
