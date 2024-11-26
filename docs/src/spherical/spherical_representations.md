# Representing Antenna Fields with Spherical Mode Expansions

One of the most popular ways to represent an antenna field is by a [spherical mode expansion](@ref sphericalexpansion).

In `AntennaFieldRepresentations.jl`, spherical expansions are represented by a struct `SphericalWaveExpansion{P, H, C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                                 | Short Description                                                |
| :---------------------------------------- | :--------------------------------------------------------------- |
| `P <: PropagationType`                    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `H <: AbstractSphericalCoefficients{C}}`  | [Type of spherical expansion coefficients](@ref abstractsphericalcoeficients). Can be [`SphericalCoefficients{C}`](@ref abstractsphericalcoeficients) for general spherical coefficients or [`FirstOrderSphericalCoefficients{C}`](@ref abstractsphericalcoeficients) for first-order spherical coefficients |
| `C <: Complex`                            | Element type of the coefficient vector                           |



!!! note
    Very efficient algorithms - e.g., for the [transmission](@ref transmission) - arise when a[`SphericalWaveExpansion`](@ref) is paired with a [`SphericalFieldSampling`](@ref).
---

## Constructors for a `SphericalWaveExpansion`
To generate a `SphericalWaveExpansion`, use one of the following constructors:

```julia
SphericalWaveExpansion{P, H, C}(coefficients::H, wavenumber <: Number) where{P <: PropagationType, C <: Complex, H<: AbstractSphericalCoefficients{C}}
```
```julia
SphericalWaveExpansion(::P, coefficients::AbstractVector{C}, wavenumber) where{P <: PropagationType, C}
```
The input arguments for the costructors are

- `P <: PropagationType` : Can be `Radiated`, `Absorbed`, or `Incident`
- `coefficients`: Coefficient vector. Can be either a struct of [`SphericalCoefficients{C}`](@ref abstractsphericalcoeficients), a struct of [`FirstOrderSphericalCoefficients{C}`](@ref abstractsphericalcoeficients), or an `AbstractVector{C}`
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## [The `AbstractSphericalCoefficients{C}` Type](@id abstractsphericalcoeficients)
As described in more detail in the [theory section](@ref sphericalexpansion), a spherical mode expansion is characterized by the expansion coefficents ``\alpha_{s \ell m}`` with ``s \in \{1,2\}``, ``\ell = 1,\ldots, L``, and ``m = -\ell, \ldots, \ell``.

In the numerical implementation, the coefficients ``\alpha_{s \ell m}`` are stored in a vector, where the triple index
``s \ell m`` is mapped to a single consecutively running index ``j`` by the function [`sℓm_to_j(s,ℓ,m)`](@ref). The inverse map from a single index to the triple index ``s \ell m`` is implemented in the function [`j_to_sℓm(j)`](@ref). 
This behaviour is implemented in implementations of the abstract type `AbstractSphericalCoefficients{C}`.
Instances of an `AbstractSphericalCoefficients{C}` can be accessed like an AbstractVector{C} with a single index `j` or with a triple index `(s,ℓ,m)`, where `j = 2 * (ℓ * (ℓ + 1) + m - 1) + s` as in the following example

Example:            
```jldoctest sphericalcoefficients
julia> using AntennaFieldRepresentations

julia> sph_coefficients= SphericalCoefficients(ComplexF64.(collect(1:16)))
16-element SphericalCoefficients{ComplexF64}:
  1.0 + 0.0im
  2.0 + 0.0im
  3.0 + 0.0im
  4.0 + 0.0im
  5.0 + 0.0im
  6.0 + 0.0im
  7.0 + 0.0im
  8.0 + 0.0im
  9.0 + 0.0im
 10.0 + 0.0im
 11.0 + 0.0im
 12.0 + 0.0im
 13.0 + 0.0im
 14.0 + 0.0im
 15.0 + 0.0im
 16.0 + 0.0im

julia> println(sph_coefficients[9]) # Here the coefficient vector is accessed with a single index
9.0 + 0.0im

julia> println(sph_coefficients[1, 2, -1]) # Here the same element is accessed with a triple index
9.0 + 0.0im

julia> sph_coefficients[1, 2, -1] = 111.0im; # setindex!() is also defined

julia> println(sph_coefficients[9])
0.0 + 111.0im
```

Besides the general `SphericalCoefficients{C}`, also the type `FirstOrderSphericalExpansion{C}` exists in `AntennaFieldRepresentations.jl`. A `FirstOrderSphericalExpansion{C}` indicates that all spherical expansion coefficients ``\alpha_{s \ell m}`` with ``|m| \neq 1`` are equal to zero. First-order spherical expansions do play a special role in the development of fast algorithms for spherical wave expansions and, thus, have their own representation in `AntennaFieldRepresentations.jl`. When a `FirstOrderSphericalExpansion{C}` is generated from an `AbstractVector{C}`, all elements which do not correspond to a triple index with ``|m|=1`` are ignored. 

Use the following constructors to generate a struct of `AbstractSphericalCoefficients{C}`
```julia
SphericalCoefficients(v::AbstractVector{C}) where{C <: Number}
```
```julia
FirstOrderSphericalCoefficients(v::AbstractVector{C}) where{C <: Number}
```

!!! note
    Since a vector of type `SphericalCoefficients` is itself a subtype of `AbstractVector`, we can easily convert a vector of `SphericalCoefficients` into a vector of first-order spherical coefficients via 
    ```julia 
    FirstorderSphericalCoefficients(x::SphericalCoefficients)
    ``` 

    Converting a an arbitrary `AbstractVector` into a `FirstOrderSphericalCoefficients`-vector drops all entries which correspond to non-first-order modes.    
---

```jldoctest sphericalcoefficients
julia> fo_sph_coefficients= FirstOrderSphericalCoefficients(sph_coefficients)
16-element FirstOrderSphericalCoefficients{ComplexF64}:
  1.0 + 0.0im
  2.0 + 0.0im
  0.0 + 0.0im
  0.0 + 0.0im
  5.0 + 0.0im
  6.0 + 0.0im
  0.0 + 0.0im
  0.0 + 0.0im
  0.0 + 111.0im
 10.0 + 0.0im
  0.0 + 0.0im
  0.0 + 0.0im
 13.0 + 0.0im
 14.0 + 0.0im
  0.0 + 0.0im
  0.0 + 0.0im
```

## Spherical Expansion Examples 

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/spherical_expansion_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/spherical_expansion_light.png" width="750" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The dipoles are visualized as arrows in this example.
  </figcaption>
</figure>
<br/>
```
