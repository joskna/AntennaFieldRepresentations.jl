# Representing Antenna Fields with Equivalent Dipole Distributions

One of the simplest ways to represent an antenna field is by a collection of electrically short (i.e., ``\ell \ll \lambda``) dipole antennas. 
Since the [radiated fields of short dipole antennas are known analytically](@ref dipole_radiated), one can simply superimpose the effects of several spatially distributed dipole antennas to approximate the radiated fields of an antenna.

In `AntennaFieldRepresentations.jl`, collections (or arrays) of electrically short dipole arrays are stored in a struct `DipoleArray{P,E,C,T}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                 | Short Description                                                |
| :------------------------ | :--------------------------------------------------------------- |
| `P <: PropagationType`    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `E <: ElmagType`          | Can be `Electric` or `Magnetic`                                  |
| `C <: Complex`            | Element type of the coefficient vector                           |
| `T <: Real`               | Number type used in the vector defining the positions of dipoles |

For extra convenience, the type aliases `HertzArray{C, T} = DipoleArray{Radiated, Electric, C, T}` and `FitzgeraldArray{C, T} = DipoleArray{Radiated, Magnetic, C, T}` are introduced. Therefore, the user will mostly interact with `HertzArray`s and `FitzgeraldArray`s while the `DipoleArray` type is hidden under the hood.

## Constructors for a `DipoleArray`
To generate a `DipoleArray`, use one of the following constructors:

```julia
DipoleArray{P, E}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
DipoleArray{P, E, C, T}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
HertzArray{C, T}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
FitzgeraldArray{C, T}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
HertzArray(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
FitzgeraldArray(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```

The input arguments for the costructors are

- `positions::Vector{V1}` : A vector of 3D-position vectors. Julia must be able to convert the type `V1` into an `SVector{3}`. Must have the same length as `orientations` and `dipolemoments`.
- `orientations::Vector{V2}`: A complex valued vector of 3D-orientations. Complex values account for elliptical polarizations in general. Julia must be able to convert the type `V2` into an `SVector{3}`. Must have the same length as `positions` and `dipolemoments`.
- `dipolemoments::Vector{C}`: A vector of complex values to denote the excitation of each individual dipole. Must have the same length as `positions` and `orientations`.
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## Dipoles with Alternative Propagation Types
Most users will probably be familiar with radiating dipoles. They correspond to the `Radiated` propagation type. 
The `Absorbed` propagation type in some sense reverses the arrow of time[^1]. Instead of radiating power away from the dipoles towards infinity, the electromagnetic fields of an `Absorbed` propagation type bring energy from infinty towards the dipole locations. 

!!! warning
    A `DipoleArray` of `Absorbed` type must not be confused with a receiving antenna!

    The `DipoleArray` is an equivalent representation of the electromagnetic fields.
    Use the [`ProbeAntenna`](@ref probeantenna) type to indicate a receiving antenna.
---

One of the main use cases of `AntennaFieldRepresentation`s of `Absorbed` type is to represent scattered fields as a superposition of `Absorbed` and `Radiated` types. The fields of `DipoleArrays` of `Absorbed` and `Radiated` type become singular at the spots where the individual dipoles are located. 

The third type of `AntennaFieldRepresentation`, i.e., the `Incident` type[^2], does not have any singularities anywhere. It can be used to represent source-free solutions of Maxwell's equations and is well suited to represent incident fields. Thus, if field representations of `Absorbed` type are not your cup of tea, you can represent any scattered field as a superposition of `Incident` and `Radiated` types.


## Dipole Examples 

!!! todo
    Coming soon ...
---

[^1]: Replacing a `Radiated` type `AntennaFieldRepresentation` by an `Absorbed` one, the electromagnetic fields of the two representations are not exactly "reversed" in time, as also the sign of the magnetic field changes. To be technically correct, both types of field representations should be considered as separate solutions of Maxwell's equations with different asymptotic boundary conditions at infinity. The fields of `DipoleArray`s of `Absorbed` type are derived from the scalar Green's function ``\mathrm{e}^{\, \mathrm{j} k r} / (4 \pi r)`` (as opposed to ``\mathrm{e}^{- \mathrm{j} k r} / (4 \pi r)`` for `Radiated` representations).

[^2]: Formally the fields of `DipoleArray`s of `Incident` type are derived from the scalar "Green's function" (more of a _pseudo_ _Green's_ _function_) ``\mathrm{sin}({\mathrm{j} k r}) / (4 \pi r)``. You can see that `Incident` fields are nothing but a superposition of `Absorbed` and `Radiated` fields.