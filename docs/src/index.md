# AntennaFieldRepresentations.jl

This package provides methods to convert between different electromagnetic field representations and to calculate the interaction/transmission between them.

!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` and SI units are used everywhere.
---
## Installation

Install `AntennaFieldRepresentations` by cloning the package to your ~\\.julia\dev folder and registering the package in the package manager (enter `]` at the julia REPL) and type 

```
pkg> add https://github.com/joskna/AntennaFieldRepresentations.jl 
```
to register the package to your `Julia` environment. Make sure that the required related packages are available.
## Featured Field Representations
The following types of electromagnetic field representations are provided by `AntennaFieldRepresentations`:

- **Dipole fields**

- **Spherical vector modes**

- **Plane wave representations**

- **Equivalent surface current**


## Operations on Field Representations
The main purpose of `AntennaFieldRepresentations` is to represent electromagnetic fields by different field expansions, manipulate the field representations, provide methods to convert different field representations into another, and to calculate the electromagnetic interactions between them. The following methods are provided by `AntennaFieldRepresentations` for physically meaningful combinations of field representations (see [Conversions Between Field Representations](@ref), [Transmission Between Field Representations](@ref), and [API Reference](@ref) for details):

-TODO: add link to convertrepresentation


-TODO: add link to convertrepresentation

Additionally, the following operations are provided (see the documentation of the individual packages for details):

- Functions to evaluate the fields at a given location
TODO: add link to field evaluations

- Functions to change into a rotated/translated coordinate system
TODO: add link to operations

- Functions to (help) transfer data into another representation with identical data structure
TODO: add link to typerelevant functions

## Related Packages
- `BEAST`
- `ClusterTrees` 

### Planned Features
- ⏳[^1] **Scalar spherical mode expansions** 
- ⏳[^1] **Inhomogeneous plane wave expansions**
- ⏳[^1] **Cylindrical vector modes**
- ⏳[^1] **Spheroidal wave functions**
- ⏳[^1] Lossless dielectrics (``\varepsilon_\mathrm{r},\mu_{\mathrm{r}} \in \mathbb{R}``)
- ⏳[^1] Lossy dielectrics (``\varepsilon_\mathrm{r},\mu_{\mathrm{r}} \in \mathbb{C}``)


## [Footnotes](@id refs)
[^1]: ⏳ Planned but not yet implemented