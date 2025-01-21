# Representing Antenna Fields with Plane-Wave Expansions

A common representation of antenna fields is in terms of [plane wave expansions](@ref planewavetheory).
In `AntennaFieldRepresentations.jl`, plane-wave expansions are represented by a struct `PlaneWaveExpansion{P, S, C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                                 | Short Description                                                |
| :---------------------------------------- | :--------------------------------------------------------------- |
| `P <: PropagationType`                    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `C <: Complex`                            | Element type of the coefficient vector                           |
| `S <: SphereSamplingStrategy, C}`         | Struct defining the sampling of the ``k``-space unit sphere. Can be [`RegularθRegularϕSampling`](@ref regularspheresampling) or [`GaussLegendreθRegularϕSampling`](@ref gausslegendresampling).|

A `PlaneWaveExpansion` of `Radiated` type corresponds to a [far field pattern ``\bm{F}(\vartheta, \varphi)``](@ref planewave_radiated) and a `PlaneWaveExpansion` of `Incident` type corresponds to a [plane wave spectrum ``\bm{P}(\vartheta, \varphi)``](@ref planewave_incident) as described in more detail in the [theory section](@ref planewavetheory).

## Constructors for a `PlaneWaveExpansion`
To generate a `PlaneWaveExpansion`, use the following constructor:
```julia
PlaneWaveExpansion(P::PropagationType, samplingstrategy::S, Eθ::Matrix{C}, Eϕ::Matrix{C}, wavenumber::Number) where{S <: SphereSamplingStrategy, C}
```
The input arguments for the constructor are

- `P <: PropagationType` : Can be `Radiated`, `Absorbed`, or `Incident`
- `samplingstrategy::SphereSamplingStrategy`: Struct defining the sampling of the ``k``-space unit sphere. Can be [`RegularθRegularϕSampling`](@ref regularspheresampling) or [`GaussLegendreθRegularϕSampling`](@ref gausslegendresampling).
- `Eθ::Matrix{C}`: θ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `Eϕ::Matrix{C}`: ϕ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## [Stored Samples of a `PlaneWaveRepresentation`](@id planewave_sampling)

In a `PlaneWaveExpansion`, an antenna far field pattern ``\bm{F}(\vartheta, \varphi)`` or a plane wave spectrum ``\bm{P}(\vartheta, \varphi)`` is representated by a discrete set of samples in the ``\vartheta, \varphi``-domain [^1] . 

The samples should be chosen such that
- integration over the sphere is fast and accurate
- local interpolation is sufficently accurate
- the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` has small computational cost (ideally should not rely on interpolation)
- redundancy in the representation is kept minimal

By default, any far-field pattern or plane-wave spectrum with a corresponding spherical mode order of ``L`` is sampled according to a [`GaussLegendreθRegularϕSampling`](@ref gausslegendresampling)
with
```math
N_\vartheta=L+1
```
and
```math
N_\varphi=2L+2 \,.
```
Although the ``\varphi``-sampling is slightly redundant (``2L+1`` samples would suffice) the even number of  ``\varphi``-samples ensures that for every sampling point there is another sampling point exactly in the oppsoite direction. This has obvious benefits for the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` as well as for pattern interpolation, when the interpolation scheme requires samples at either side of the coordinate poles at ``\vartheta\in \{0,\pi\}``.

!!! tip
    `PlaneWaveExpansion`s can be sampled according to any arbitrary [`SphereSamplingStrategy`](@ref spheresampling).
 ---


## Methods Special to `PlaneWaveExpansion`s
In addition to the methods defined in the interface of `AntennaFieldRepresentation`, a `PlaneWaveExpansion` supports the following methods:

!!! todo
    Add a comprehensive list of `PlaneWaveExpansion` methods
---

## [References](@id planewave_refs)
[^1]: Since the propagation vector ``\bm{k}`` and the corresponding unit vector into the same direction ``\hat{\bm{k}}`` encode the propagation direction of the plane wave, we may use the shorthand notation ``\bm{P}(\hat{\bm{k}})`` to represent the slightly longer expression ``\bm{P}(\vartheta, \varphi)`` whenever convenient (sometimes we may even mix ``\vartheta, \varphi`` with ``\hat{\bm{k}}`` in the same expression). 
This should not lead to any ambiguities because the relation between ``\hat{\bm{k}}`` and the tuple ``\vartheta, \varphi`` is one-to one, as each encodings uniquely define the same point on the unit sphere.