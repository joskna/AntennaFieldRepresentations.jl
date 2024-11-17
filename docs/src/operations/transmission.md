# [Transmission between an `AntennaFieldRepresentation` and a `FieldSampling`](@id transmission) 
The [`transmit`](@ref) method can be used to evaluate the probe outputs of a [`FieldSampling`](@ref fieldsampling) for a given [`AntennaFieldRepresentation`](@ref) without constructing a [`TransmitMap`](@ref) first. However, using the `transmit` method is only recommended if the operation is not called multiple times. In most cases, it is advisable to construct and use a `TransmitMap` to avoid redundant allocations. 

## [`TransmitMap`s](@id transmissionmap)
A `TransmitMap` is constructed via the constructor

```julia
A = TransmitMap(aut_field::AntennaFieldRepresentation, fieldsampling::Fieldsampling)
```

Here is an example of how to construct a `TransmitMap` between a `SphericalWaveExpansion` and a `SphericalFieldSampling`:
```jldoctest transmitmaps; output=false
using AntennaFieldRepresentations

Z₀=376.730313669;                       # free-space wave impedance
f= 1.5e9;                               # frequency
λ = AntennaFieldRepresentations.c₀ / f; # wavelength
k0 = 2 * pi / λ;                        # wavenumber

sph_coefficients = SphericalCoefficients(ComplexF64.(collect(1:16)))
aut_field = SphericalWaveExpansion(Radiated(),sph_coefficients, k0)

Lmax = 2 # maximum mode order is 2 in this example

# spherical coefficients of incident field Hertzian dipole at 1.0m distance
αinc=AntennaFieldRepresentations.αinc_dipole(1.0, Lmax, k0)

Jθ= 2*Lmax + 1;
Jϕ= 2*Lmax + 1;
samplingstrategy= RegularθRegularϕSampling(Jθ, Jϕ) 

fieldsampling = SphericalFieldSampling(samplingstrategy, αinc);

A=TransmitMap(aut_field, fieldsampling)

# output

30×16 SphericalTransmitMap{ComplexF64}
```

Depending on the concrete types of the `aut_field` and the `fieldsampling`, different types of `TransmitMap` will be created to implement the most efficient agorithm for the transmission. The following cross table clarifies which type of `TransmitMap` is created by which pair of input types.

|                            |    | `IrregularFieldSampling` | `SphericalFieldSampling`            | 
| :------------------------  |:-- |:------------------------ | :---------------------------------- |
|                            |    |                          |                                     | 
|**`DipoleArray`**           |    | `MLFMMTransmitMap`       | `MLFMMSphericalTransmitMap`         |
|**`SphericalWaveExpansion`**|    | `DirectTransmitMap`      |[`SphericalTransmitMap`](@ref)[^1]   |
|**`PlaneWaveExpansion`**    |    | `DirectTransmitMap`      | `PlaneWaveSphericalTransmitMap`[^1] |
|**`SurfaceCurrents`**       |    | `MLFMMTransmitMap`       | `MLFMMSphericalTransmitMap`         |


[^1]: A direct (non-iterative) inverse transmit map is available using the `inverse()` command.

As explained in the [`OperationMaps` Chapter of this documentation](@ref operationmaps_linmap), each `TransmitMap` operates as linear map and its `adjoint`, `transpose` and `inverse` map can be constructed via

```jldoctest transmitmaps
julia>  Aᴴ = adjoint(A)
16×30 LinearMaps.AdjointMap{ComplexF64} of
  30×16 SphericalTransmitMap{ComplexF64}

julia>  Aᵀ = transpose(A)
16×30 LinearMaps.TransposeMap{ComplexF64} of
  30×16 SphericalTransmitMap{ComplexF64}

julia> A⁻¹ = inverse(A)
16×30 AntennaFieldRepresentations.InverseSphericalTransmitMap{ComplexF64}
```

!!! note
    Not all types of `TransmitMap` allow to setup an inverse efficiently. In most cases, a direct inverse is not available but the inverse map is implemented via an iterative solver.  The types of `TransmitMap` which have a direct inverse available, are annotated with a superscript [^1] in the table above.  
---

Since the different types of `TransmitMap` represent different algorithms for the computation of the transmission, they have different parameters which may be tuned by the user. The parameters of the `TransmitMaps` can be changed after an instance of the `TransmitMap` has been created using the command
```julia
changeparameters(A::TransmitMap; kwargs...)
```
where the available keyword arguments depend on the concrete type of the `TransmitMap`

Each concrete type of `TransmitMap` is explained in more detail. The sections briefly describe the underlying algorithms and specify the available keyword arguments to tune the parameters of the underlying algorithms.

## `SphericalTransmitMap`
The `SphericalTransmitMap` type calculates the transmission between a `SphericalWaveExpansion` and a `SphericalFieldSampling`. It implements the spherical Wacker algorithm which is described in more detail [in the electromagnetic theory section](@ref fastsphericalS12).

The keyword arguments for `changeparameters(A::SphericalTransmitMap; kwargs...)` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`firstorder::Bool` |If true, spherical coefficients of the incident probe field are treated as `FirstOrderSphericalCoefficients`. A change from `false` to `true` drops all entries which correspond to non-first-order modes. A change from `true` to `false` replaces non-first-order modes by zeros.|
|`samplingstrategy::SphereSamplingStrategy`|Replace the `SphereSamplingStrategy`. This effectively creates a new instance of `SphericalTransmitMap` from scratch.|
|`incidentcoefficients::AbstractVector`| Replace the spherical coefficients of the incident probe field by the entries of the input vector. If the `SphericalTransmitMap` is `firstorder`, all entries of the input vector which correspond to non-first-order modes are ignored.|

## `PlaneWaveSphericalTransmitMap`
The `PlaneWaveSphericalTransmitMap` type calculates the transmission between a `PlaneWaveExpansion` and a `SphericalFieldSampling`.
It is implemented as a concatination of a `PlaneWaveToSphericalMap` and a `SphericalTransmitMap`, i.e., the `PlanewaveExpansion` is converted into a `SphericalWaveExpansion` and the transmission is executed as a `SphericalTransmitMap`.

The keyword arguments for `changeparameters(A::PlaneWaveSphericalTransmitMap; kwargs...)` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`firstorder::Bool` |If true, spherical coefficients of the incident probe field are treated as `FirstOrderSphericalCoefficients`. A change from `false` to `true` drops all entries which correspond to non-first-order modes. A change from `true` to `false` replaces non-first-order modes by zeros.|
|`samplingstrategy::SphereSamplingStrategy`|Replace the `SphereSamplingStrategy`. This effectively creates a new instance of `SphericalTransmitMap` from scratch.|
|`incidentcoefficients::AbstractVector`| Replace the spherical coefficients of the incident probe field by the entries of the input vector. If the `SphericalTransmitMap` is `firstorder`, all entries of the input vector which correspond to non-first-order modes are ignored.|

## `MLFMMTransmitMap`

## `MLFMMSphericalTransmitMap`
The `MLFMMSphericalTransmitMap`type calculates the transmission between an `AntennaFieldRepresentation` based on equivalent currents and a `SphericalFieldSampling`.
In its implementation, the far-field pattern (i.e., a `PlaneWaveExpansion`) is calculated from the equivalent current representation first. For maximum efficiencyc the source region is hierrchically divided into boxes arranged in an octree. The overall far-field pattern is calculated by aggregating the far-fields of sources in smaller boxes in a [MLFMM-scheme](@ref mlfmm_aggregation).
In the second step, the transmission between the calculated `PlaneWaveExpansion` and the `Spherical FieldSampling` is calculated through a `PlaneWaveSphericalTransmitMap`.

The keyword arguments for `changeparameters(A::MLFMMSphericalTransmitMap; kwargs...)` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`firstorder::Bool` |If true, spherical coefficients of the incident probe field are treated as `FirstOrderSphericalCoefficients`. A change from `false` to `true` drops all entries which correspond to non-first-order modes. A change from `true` to `false` replaces non-first-order modes by zeros.|
|`samplingstrategy::SphereSamplingStrategy`|Replace the `SphereSamplingStrategy`. This effectively creates a new instance of `SphericalTransmitMap` from scratch.|
|`incidentcoefficients::AbstractVector`| Replace the spherical coefficients of the incident probe field by the entries of the input vector. If the `SphericalTransmitMap` is `firstorder`, all entries of the input vector which correspond to non-first-order modes are ignored.|

## `DirectTransmitMap`

