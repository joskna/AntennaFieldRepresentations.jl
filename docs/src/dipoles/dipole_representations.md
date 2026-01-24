# Representing Antenna Fields with Equivalent Dipole Distributions

One of the simplest ways to represent an antenna field is by a collection of electrically short (i.e., ``\ell \ll \lambda``) dipole antennas. 
Since the [radiated fields of short dipole antennas are known analytically](@ref dipole_radiated), one can simply superimpose the effects of several spatially distributed dipole antennas to approximate the radiated fields of an antenna.

In `AntennaFieldRepresentations.jl`, collections (or arrays) of electrically short dipole arrays are stored in a struct `DipoleArray{P,E,T,C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                 | Short Description                                                |
| :------------------------ | :--------------------------------------------------------------- |
| `P <: PropagationType`    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `E <: ElmagType`          | Can be `Electric` or `Magnetic`                                  |
| `T <: Real`               | Number type used in the vector defining the positions of dipoles |
| `C <: Complex`            | Element type of the coefficient vector                           |


For extra convenience, the type aliases `HertzArray{T,C} = DipoleArray{Radiated, Electric, T, C}` and `FitzgeraldArray{T, C} = DipoleArray{Radiated, Magnetic, T, C}` are introduced. Therefore, the user will mostly interact with `HertzArray`s and `FitzgeraldArray`s while the `DipoleArray` type is hidden under the hood.

## Constructors for a `DipoleArray`
To generate a `DipoleArray`, use one of the following constructors:

```julia
DipoleArray{P, E}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
DipoleArray{P, E, T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
HertzArray{T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
FitzgeraldArray{T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
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

Let us first create an Array of Hertzian dipoles. 
```jldoctest dipoleexamples ; output=false
using AntennaFieldRepresentations

f = 1.5e9;  # Set frequency to 1.5 GHz
λ = AntennaFieldRepresentations.c₀ / f;  # wavelength
k0 = 2 * pi / λ;  # wavenumber

positions= [[-λ, 0, 0], [0, λ/2, λ/2], [0, -λ,0]];
orientations= [complex.([0.0,0.0,1.0]), complex.([0.0,1.0,0.0]), complex.([0.1,0.0,0.0])];
dipolemoments= [ComplexF64(1.0), ComplexF64(1.0), ComplexF64(1.0)];

dipoles = HertzArray(positions, orientations, dipolemoments, k0);

# output
3-element HertzArray{Float64, ComplexF64}:
 1.0 + 0.0im
 1.0 + 0.0im
 1.0 + 0.0im

```

The resulting set of dipoles might be visualized as follows:

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/dipoles_light.png" width="750" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The dipoles are visualized as arrows in this example.
  </figcaption>
</figure>
<br/>
```

We might be interested in the electromagnetic field - maybe its ``E_x`` component - which is radiated by these dipoles, let's say in a plane ``z=5\lambda``.
Thus, we can define an array of observation points and evaluate the electric field at these observation points as follows (the `Ref()` command ensures that this input is treated as a constant for Julia's broadcasting operator "`.`"):
```jldoctest dipoleexamples ; output=false
Rs= [[x , y, 5λ] for x in -10λ:λ/4:10λ , y in -10λ:λ/4:10λ] # Define observation points

E = efield.(Ref(dipoles), Rs) # Evaluate E-field at observation points
Ex=[e[1] for e in E] # Extract x-component of E-field
size(Ex)

# output

(81, 81)

```

The resulting field (stored as an ordinary matrix) can then be visualized, e.g., with `Makie.jl` or `Plots.jl`
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/field_dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/field_dipoles_light.png" width="750" >
  <img alt="" src="" width="750">
</picture>

  <figcaption>
    Magnitude of the Eₓ-component of the radiated field in a plane at z=5λ.
  </figcaption>
</figure>
<br/>
```

Furthermore, we might be interested in the far fields radiated by the dipole collection. 
Thus, we define pairs `(θ, ϕ)` of angles for the directions in which we want to evaluate the far fields and calculate the far field via the `farfield` command 
```jldoctest dipoleexamples ; output=false
θs= LinRange(0.0, π, 40)
ϕs= LinRange(0.0, 2π, 80)

directions= [(θ, ϕ) for θ in θs, ϕ in ϕs]

farfields= farfield.(Ref(dipoles), directions)
size(farfields)

# output

(40, 80)

```
The output is a pair `(Fθ, Fϕ)` for each input in the `directions` array. We can retrieve the individual components via
```jldoctest dipoleexamples ; output=false
Fθ= [ff[1] for ff in farfields]
Fϕ= [ff[2] for ff in farfields]
size(Fϕ)

# output

(40, 80)

```

The resulting  far field can then be visualized, e.g., with `Makie.jl` or `Plots.jl`
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/farfield_dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/farfield_dipoles_light.png" width="750" >
  <img alt="" src="" width="750">
</picture>

  <figcaption>
    Magnitude of  of the radiated far field in dB-scale.
  </figcaption>
</figure>
<br/>
```

As it turns out, the radiated far-field is rather omni-directional. With deeper thought, this might not be too surprising because each field component is equally well excited by the three dipoles in ``x-``, ``y-``, and ``z-`` direction.


[^1]: Replacing a `Radiated` type `AntennaFieldRepresentation` by an `Absorbed` one, the electromagnetic fields of the two representations are not exactly "reversed" in time, as also the sign of the magnetic field changes. To be technically correct, both types of field representations should be considered as separate solutions of Maxwell's equations with different asymptotic boundary conditions at infinity. The fields of `DipoleArray`s of `Absorbed` type are derived from the scalar Green's function ``\mathrm{e}^{\, \mathrm{j} k r} / (4 \pi r)`` (as opposed to ``\mathrm{e}^{- \mathrm{j} k r} / (4 \pi r)`` for `Radiated` representations).

[^2]: Formally the fields of `DipoleArray`s of `Incident` type are derived from the scalar "Green's function" (more of a _pseudo_ _Green's_ _function_) ``\mathrm{sin}({\mathrm{j} k r}) / (4 \pi r)``. You can see that `Incident` fields are nothing but a superposition of `Absorbed` and `Radiated` fields.