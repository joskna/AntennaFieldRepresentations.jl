# [MLFMM Representations of Electromagnetic Fields](@id mlfmm_source)
 The source volume of an equivalent source representation of a certain antena field (e.g., equivalent surface currents or an array of dipoles) may extend over a large volume (i.e., the diameter exceeds several wavelength) for certain fields. This may lead to a large number of source coefficients representing the electromagnetic field.
 
 In these cases, it may be beneficial to represent the field contributions in a MLFMM datastructure which hierarchically divides the source volume into nested boxes. The `MLFMMSource` data structure implements this hierarchically nested representation scheme.

 In `AntenneFieldrepresentations.jl`, this hierarchically nested data structure is implemented by the struct `MLFMMSource{M,Y,C,X}` which is a sub type of `AntennaFieldRepresentation{Radiated,C}`. 

The type parameters have the following meaning:

| Parameter                                 | Short Description                                                |
| :---------------------------------------- | :--------------------------------------------------------------- |
| `M <: ResampleMap`                        | Type of the `ResampleMap`s which are used to interpolate from one box level to the next |
| `Y <: SphereSamplingStrategy`             | Type of the sampling strategy of the `PlaneWaveRepresentation`s used to represent the accumulated field contributions of each box |
| `C <: Complex`                            | Element type of the coefficient vector                           |
| `X <: MLFMMTree`                          | Type of the tree structure used to hirarchically split the source domain into subdomains|

!!! tip
    Field representations by a `MLFMMSource` very efficiently implement 
    - the conversion into a `PlaneWaveRepresentation{Radiated}` (i.e., a complete far-field patern)
    - the conversion into a `SphericalWaveExpansion` (via the detour over a conversion into a `PlaneWaveRepresentation{Radiated}`)
    - the transmission operator with an `IrregularFieldSampling` (as part of an `MLFMMTransmitMap`)
    - the transmission operator with a `SphericalFieldSampling` (via the detour over a conversion into a `SphericalWaveExpansion`)  
---

## Constructors for an `MLFMMSource`
To generate an `MLFMMSource`, use the following constructor:
```julia
MLFMMSource(
    basisfunctions, #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions) 
    wavenumber::T;
    expectedaccuracy = T(1e-3),
    verbose = false,
    minboxlength = π / (2*wavenumber),
    orderθ = 8,
    orderϕ = 8,
    samplingtype::Type{S} = GaussLegendreθRegularϕSampling,
) where {T<:Real,S<:SphereSamplingStrategy}
```

The input arguments for the constructor are

- `basisfunctions` : Original antenna field representation which is the basis for the `MLFMMsource`. Can be of type `SurfaceCurrentDensity{Radiated}`, `DipoleArray{Radiated}` or (will soon be implemented) `NamedTuple(:points,:sourcefunctions)`
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

The optional keyword arguments are
- `expectedaccuracy = T(1e-3)`: Defines the desired approximation accuracy for all operations involving the `MLFMMsource`. Mainly affects the number of samples used to store the `PlaneWaveRepresentation`s for each box.
- `verbose = false`: Determines if (debug) text should be printed to console for all operations involving the `MLFMMSource`.
- `minboxlength = π / (2 * wavenumber)`: Side length of the smallest box, defaults to half a wavelength. All larger box side lengths are a ``2^n``-multiple of the `minboxlength`.
- `orderθ=8`: Interpolation order along ``\theta`` for the `ResampleMap`s used to intrpolate from one box level to the next. Only affects `ResampleMaps` with local interpolation strategies. 
- `orderϕ=8`: Interpolation order along ``\phi`` for the `ResampleMap`s used to intrpolate from one box level to the next. Only affects `ResampleMaps` with local interpolation strategies.
- `samplingtype::Type{SphereSamplingStrategy} = GaussLegendreθRegularϕSampling`: Type of the sampling strategy of the `PlaneWaveRepresentation`s used to represent the accumulated field contributions of each box
