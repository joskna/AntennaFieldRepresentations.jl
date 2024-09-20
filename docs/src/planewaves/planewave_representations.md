# Representing Antenna Fields with Plane-Wave Expansions

A common representation of antenna fields is in terms of [plane wave expansions](@ref planewavetheory).
In `AntennaFieldRepresentations.jl`, plane-wave expansions are represented by a struct `PlaneWaveExpansion{P, S <: SphereSamplingStrategy, C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                                 | Short Description                                                |
| :---------------------------------------- | :--------------------------------------------------------------- |
| `P <: PropagationType`                    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `C <: Complex`                            | Element type of the coefficient vector                           |
| `S <: SphereSamplingStrategy, C}`         | Struct defining the sampling of the ``k``-space unit sphere. Can be [`RegularSphereSamplingStrategy`](@ref regularspheresampling) or [`GaussLegendreθRegularϕ`](@ref gausslegendresampling).|

## Constructors for a `PlaneWaveExpansion`
To generate a `PlaneWaveExpansion`, use the following constructor:
```julia
PlaneWaveExpansion(P::PropagationType, samplingstrategy::S, Eθ::Matrix{C}, Eϕ::Matrix{C}, wavenumber::Number) where{S <: SphereSamplingStrategy, C}
```
The input arguments for the costructor are

- `P <: PropagationType` : Can be `Radiated`, `Absorbed`, or `Incident`
- `samplingstrategy::SphereSamplingStrategy`: Struct defining the sampling of the ``k``-space unit sphere. Can be [`RegularSphereSamplingStrategy`](@ref regularspheresampling) or [`GaussLegendreθRegularϕ`](@ref gausslegendresampling).
- `Eθ::Matrix{C}`: θ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `Eϕ::Matrix{C}`: ϕ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f`

## [The `SphereSamplingStrategy` Type](@id planewave_sampling)

In a `PlaneWaveExpansion`, an antenna far field pattern ``\bm{F}(\vartheta, \varphi)`` or a plane wave spectrum ``\bm{P}(\vartheta, \varphi)`` is representated by a discrete set of samples in the ``\vartheta, \varphi``-domain. 

The samples should be chosen such that
- integration over the sphere is fast and accurate
- local interpolation is sufficently accurate
- the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` has small computational cost (ideally should not rely on interpolation)
- redundancy in the representation is kept minimal

### [`RegularSphereSamplingStrategy`](@id regularspheresampling)

### [`GaussLegendreθRegularϕ`](@id gausslegendresampling)
Expansions with a `GaussLegendreθRegularϕ` sampling are sampled along ``\varphi`` with  ``2L+2`` equally distributed samples  in ``0\leq \varphi < 2\pi`` and distribute the ``\vartheta``-samples on a Gauß-Legendre based grid with ``0< \vartheta < \pi``, according to 
```math
\vartheta_k=\text{arccos}(-x_k)
```
where ``x_k`` are the roots of the ``n``th Legendre polynomial. To represent any far-field pattern or plane-wave spectrum up to a spherical mode order of ``L``, the number of ``\vartheta``- and ``\varphi``-samplesis 
```math
N_\vartheta=L+1
```
```math
N_\varphi=2L+2
```
Although the ``\varphi``-sampling is slightly redundant (``2L+1`` samples would suffice) the even number of  ``\varphi``-samples ensures that for every sampling point there is another sampling point exactly in the oppsoite direction. This has obvious benefits for the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` as well as for pattern interpolation, when the interpolation scheme requires samples at either side of the coordinate poles at ``\vartheta\in \{0,\pi\}``.
The figures below compare the regular sampling points in a `RegularSphereSamplingStrategy` and the Lagrange sampling points utilized by `GaussLegendreθRegularϕ` sampling. 


```@raw html
<figure>
<picture>
  <source  srcset="../assets/sampling_regular.svg" width="200">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    This is a regularly sampled sphere.
  </figcaption>
</figure>
<br/>
```
```@raw html
<figure>
<picture>
  <source  srcset="../assets/sampling_lagrange.svg" width="200">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    This sphere is sampled according to a Gauß-Legendre sampling rule.
  </figcaption>
</figure>
<br/>
```

The sampling choice naturally leads to an accurate integration rule over the complete Ewald sphere 
```math
\oiint f(\hat{\bm{k}}) \, \mathrm{d}^2 \hat{\bm{k}} 
= 
\int \limits_{0}^{2\pi} \int \limits_{0}^{\pi} \, f(\vartheta, \varphi) \, \sin{\vartheta}\, \mathrm{d}\vartheta
\mathrm{d} \varphi 
\approx 
\sum \limits_{i=1}^{N_\varphi} \sum \limits_{k=1}^{N_\vartheta} \dfrac{2\pi}{N_\varphi}
w_k\, f(\vartheta_k, \varphi_i)
``` 
where ``w_k`` are the Gauß-Legendre quadrature weights. The calculation of the Gauß-Legendre quadrature weights ``w_k`` and quadrature nodes ``x_k`` is efficiently implemented (to compute the ``N_\vartheta``-point Gauß quadrature nodes and weights to 16-digit accuracy takes ``\mathcal{O}(N_\vartheta)`` time [^1]) in the Julia package [`FastGaussQuadrature.jl`](https://juliaapproximation.github.io/FastGaussQuadrature.jl/stable/).


## Methods Special to `PlaneWaveExpansion`s
In addition to the methods defined in the interface of `AntennaFieldRepresentation`, a `PlaneWaveExpansion` supports the following methods:

## [References](@id planewave_refs)
[^1]: I. Bogaert, Iteration-free computation of Gauss–Legendre quadrature nodes and weights, SIAM J.ournal on Scientific Computing, Vol. 36, No. 3, pp. A1008–A1026, 2014 [DOI. 10.1137/140954969](https://doi.org/10.1137/140954969)