<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.svg" height="190">
  <img alt="" src="" height="190">
</picture>

[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joskna.github.io/AntennaFieldRepresentations.jl/dev/)

# AntennaFieldTransformations

## Introduction
This package contains various representations of electromagnetic fields and their interactions with antennas. 
The goal is to provide (efficient) methods to transfrom any representation into another to be able to choose the most suitable representation for the problem at hand.
The purpose of this package is primarily academic and educational. At the moment and the package is under heavy development and represents by no means a finished state. Use with care!

## Installation
The package is not registered in the Julia General Registry, yet, but a future version of this package might.
For the time being, you can install the package by switching into the package manager (type `]` in the REPL) and run
```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl
```
