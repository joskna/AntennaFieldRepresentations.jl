# Representation of Antenna Fields

The electromagnetic fields of an antenna can be represented in various ways. Some of the most important field representations are implemented in the package `AntennaFieldRepresentations.jl`. 

The implemented field representations are `struct`s which are subtypes of the abstract supertype `AntennaFieldRepresentation`. 
An `AntennaFieldRepresentation` is interpreted as a vector of coefficients plus additional context which allows to compute the electromagnetic fields from the coefficient vector. As a consequence of the interpretation of all `AntennaFieldrepresentation`s as coefficient vectors, the `AntennaFieldRepresentation`-type is a subtype of Julia's `Base.AbstractVector`.

For all `AntennaFieldRepresentation`s, one may evaluate their E-field or H-field (or both) at a location ``\bm{R}`` in space with one of the functions `efield`, `hfield`, or `ehfield`, respectively.

!!! tip
    All implementations of the type `AntennaFieldrepresentation` support the interface of `Base.AbstractVector`. In particular the methods `Base.getindex` and `Base.setindex!` are supported which allow to access the coefficient vector by indexing the struct representing the fields directly by square brackets `[]`. Furthermore, the functions `Base.size` and `Base.similar` are supported, as well.
---

