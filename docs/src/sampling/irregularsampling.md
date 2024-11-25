# [Irregularly Distributed Field Samplings](@id irregularsampling)
The most general [`FieldSampling`](@ref fieldsampling) in `AntennaFieldRepresentations.jl` is the `IrregularFieldSampling`. 
It allows to sample the antenna field at arbitrary positions with arbitrary probe antennas.
The user must specify
- a list of measurement positions 
- a vector of [probe antennas](@ref probeantenna) which are used in the measurement scenario
- a list of probe IDs to assign the utilized probe to each measurement position
- a list of Euler angles ``\vartheta``, ``\varphi``, ``\chi`` to define the probe rotation at each measurement position

## [Probe Antenna](@id probeantenna)
For the definition of a probe antenna, `AntennaFieldRepresentations.jl` provides the type [`ProbeAntenna`](@ref). The `ProbeAntenna` type is a thin wrapper around any radiating [`AntennaFieldRepresentation`](@ref fieldrepresentation). A `ProbeAntenna` object can be generated by `ProbeAntenna(aut_field, size)` where 
- `aut_field <: AntennaFieldRepresentation{Radiating, C<: ComplexF64}` may be any[`AntennaFieldRepresentation`](@ref fieldrepresentation) of radiating type 
- `size <: Real` is a real number denotig the physical size of the probe antenna in terms of the radius of the smallest sphere enclosing the complete probe antenna.

## Constructors of `IrregularFieldSampling`
!!! note
    Coming soon ...
---

## E-Field and H-Field Requests
Since unfiltered sampling of the electric or magnetic field of an [`AntennaFieldRepresentation`](@ref fieldrepresentation) at given measurement locations is a common use case, a convenience constructor is provided by `AntennaFieldRepresentations.jl` for these cases. The methods `EfieldSampling(positions::Vector{V})` and `HfieldSampling(positions::Vector{V})` return `IrregularFieldSampling`s which correspond to sampling the field of an `AntennaFieldRepresentation` at given positions with [Hertzian dipole probes](@ref receivedsignaldipoles) or [Fitzgerald dipole probes](@ref receivedsignaldipoles) respectively.