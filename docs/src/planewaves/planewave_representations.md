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
The input arguments for the costructor are

- `P <: PropagationType` : Can be `Radiated`, `Absorbed`, or `Incident`
- `samplingstrategy::SphereSamplingStrategy`: Struct defining the sampling of the ``k``-space unit sphere. Can be [`RegularθRegularϕSampling`](@ref regularspheresampling) or [`GaussLegendreθRegularϕSampling`](@ref gausslegendresampling).
- `Eθ::Matrix{C}`: θ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `Eϕ::Matrix{C}`: ϕ-component amplitudes of the plane waves. Dimensions must match the `samplingstrategy`.
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## [The `SphereSamplingStrategy` Type](@id planewave_sampling)

In a `PlaneWaveExpansion`, an antenna far field pattern ``\bm{F}(\vartheta, \varphi)`` or a plane wave spectrum ``\bm{P}(\vartheta, \varphi)`` is representated by a discrete set of samples in the ``\vartheta, \varphi``-domain [^1] . 

The samples should be chosen such that
- integration over the sphere is fast and accurate
- local interpolation is sufficently accurate
- the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` has small computational cost (ideally should not rely on interpolation)
- redundancy in the representation is kept minimal

### [`RegularθRegularϕSampling`](@id regularspheresampling)
Objects with a `RegularθRegularϕSampling` are sampled along ``\varphi`` with ``N_\varphi`` equally distributed samples  in ``0\leq \varphi < 2\pi`` and along ``vartheta`` with ``N_\vartheta`` equally distributed samples  in ``0\leq \varphi \leq 2\pi``.

The sampling steps ``\Delta\vartheta = 2\pi / J_\vartheta`` and ``\Delta\varphi = 2\pi / J_\varphi`` must be whole-number fractions of ``2\pi`` to guarantee that the sampling step between the two samples at the end of a ``\varphi``-ring or a ``\vartheta``
ring is the same as everywhere else.
This is enforced during the construction of a `RegularθRegularϕSampling` as the whole-number divisors ``J_\vartheta`` and ``J_\varphi`` are the input arguments for the constructor.

To generate a `RegularθRegularϕSampling` struct, use the following constructor
```julia
RegularθRegularϕSampling(Jθ::Integer, Jϕ::Integer)
```

While local and even global interpolation is easy and efficient with a `RegularθRegularϕSampling`, the drawback is that integration is a bit more costly as compared to a `GaussLegendreθRegularϕSampling`.

### [`GaussLegendreθRegularϕSampling`](@id gausslegendresampling)
Objects with a `GaussLegendreθRegularϕSampling` are sampled along ``\varphi`` with ``N_\varphi`` equally distributed samples  in ``0\leq \varphi < 2\pi`` and ``N_\vartheta`` ``\vartheta``-samples on a Gauß-Legendre based grid with ``0< \vartheta < \pi``, according to 
```math
\vartheta_k=\text{arccos}(-x_k)
```
where ``x_k`` are the roots of the ``N_\vartheta``th Legendre polynomial. 

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
where ``w_k`` are the Gauß-Legendre quadrature weights. The calculation of the Gauß-Legendre quadrature weights ``w_k`` and quadrature nodes ``x_k`` is efficiently implemented (to compute the ``N_\vartheta``-point Gauß quadrature nodes and weights to 16-digit accuracy takes ``\mathcal{O}(N_\vartheta)`` time [^2]) in the Julia package [`FastGaussQuadrature.jl`](https://juliaapproximation.github.io/FastGaussQuadrature.jl/stable/).

To generate a `GaussLegendreθRegularϕSampling` struct, use the following constructor
```julia
GaussLegendreθRegularϕSampling(Nθ::Integer, Nϕ::Integer)
```
where `Nθ` corresponds to the stored number of samples in ``\vartheta``-direction and `Nϕ` corresponds to the number of samples in ``\varphi``
direction.

To represent any far-field pattern or plane-wave spectrum up to a spherical mode order of ``L``, the number of ``\vartheta``- and ``\varphi``-samples should be chosen according to 
```math
N_\vartheta=L+1
```
```math
N_\varphi=2L+2
```
Although the ``\varphi``-sampling is slightly redundant (``2L+1`` samples would suffice) the even number of  ``\varphi``-samples ensures that for every sampling point there is another sampling point exactly in the oppsoite direction. This has obvious benefits for the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` as well as for pattern interpolation, when the interpolation scheme requires samples at either side of the coordinate poles at ``\vartheta\in \{0,\pi\}``.
The figures below compare the regular sampling points in a `RegularθRegularϕSampling` and the Lagrange sampling points utilized by a `GaussLegendreθRegularϕSampling`. 


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

## Methods Special to `PlaneWaveExpansion`s
In addition to the methods defined in the interface of `AntennaFieldRepresentation`, a `PlaneWaveExpansion` supports the following methods:

!!! todo
    Add a comprehensive list of `PlaneWaveExpansion` methods
---

## [References](@id planewave_refs)
[^1]: Since the propagation vector ``\bm{k}`` and the corresponding unit vector into the same direction ``\hat{\bm{k}}`` encode the propagation direction of the plane wave, we may use the shorthand notation ``\bm{P}(\hat{\bm{k}})`` to represent the slightly longer expression ``\bm{P}(\vartheta, \varphi)`` whenever convenient (sometimes we may even mix ``\vartheta, \varphi`` with ``\hat{\bm{k}}`` in the same expression). 
This should not lead to any ambiguities because the relation between ``\hat{\bm{k}}`` and the tuple ``\vartheta, \varphi`` is one-to one, as each encodings uniquely define the same point on the unit sphere.

[^2]: I. Bogaert, Iteration-free computation of Gauss–Legendre quadrature nodes and weights, SIAM J.ournal on Scientific Computing, Vol. 36, No. 3, pp. A1008–A1026, 2014 [DOI. 10.1137/140954969](https://doi.org/10.1137/140954969)