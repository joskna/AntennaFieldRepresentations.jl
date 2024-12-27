# [Spherically Sampled Data Structures](@id spheresampling)

In `AntennaFieldRepresentations.jl`, several data structures represent functions ``\bm{f}(\vartheta, \varphi)`` which are sampled on a sphere at certain sampling points ``(\vartheta_i, \varphi_i)``. Examples for spherically sampled data structures are [`PlaneWaveExpansion`](@ref)s or [`SphericalFieldSampling`](@ref)s. 

A natural choice for the sampling points ``(\vartheta_i, \varphi_i)`` is to choose a tensor product of 1D samplings along ``\vartheta`` and ``\varphi``, resepectively, i.e., by defining ``N_\vartheta`` sampling points ``\vartheta_m`` along ``\vartheta`` and ``N_\varphi`` sampling points ``\varphi_n`` along ``\varphi``, we end up with ``N_\vartheta \times N_\varphi`` pairs ``(\vartheta_m, \varphi_n)`` which define the sampling points on the sphere.

Internally, the data can be stored in a matrix[^1], where the ``m``th row corresponds the angle ``\vartheta_m`` and the ``n``th coloumn corresponds to the angle ``\varphi_n``. The ``m, n``- entry of the matrix therefore corresponds to the function value ``\bm{f}(\vartheta_m, \varphi_n)``.

Different choices for the distribution of sampling locations may have application specific benefits or drawbacks. To enable the users to use the most suitable choice for their problem at hand, several different sampling distributions are implemented in `AntennaFieldRepresentations.jl`.
The values of the sampling angles are defined by a [`SphereSamplingStrategy`](@ref) object. 

Each spherically sampled data structure contains a `SphereSamplingStrategy` object such that we know which sampling locations ``(\vartheta_m, \varphi_n)`` correspond to the matrix entries.

The following `SphereSamplingStrategy`s are implemented in `AntennaFieldRepresentations.jl`:

## [`RegularθRegularϕSampling`](@id regularspheresampling)
Objects with a `RegularθRegularϕSampling` are sampled along ``\varphi`` with ``N_\varphi`` equally distributed samples  in ``0\leq \varphi < 2\pi`` and along ``\vartheta`` with ``N_\vartheta`` equally distributed samples  in ``0\leq \varphi \leq \pi``.

We have 
- ``\vartheta_m = (m-1)\, \Delta\vartheta`` with ``m=1,\ldots,N_\vartheta`` 
- ``\varphi_n = (n-1)\, \Delta\varphi`` with ``n=1,\ldots,N_\varphi``

The sampling steps ``\Delta\vartheta = 2\pi / J_\vartheta`` and ``\Delta\varphi = 2\pi / J_\varphi`` must be whole-number fractions of ``2\pi`` to guarantee that the sampling step between the two samples at the end of a ``\varphi``-ring or a ``\vartheta``
ring is the same as everywhere else.
This is enforced during the construction of a `RegularθRegularϕSampling` as the whole-number divisors ``J_\vartheta`` and ``J_\varphi`` are the input arguments for the constructor.

To generate a `RegularθRegularϕSampling` struct, use the following constructor
```julia
RegularθRegularϕSampling(Jθ::Integer, Jϕ::Integer)
```

!!! note
    You can define a `RegularθRegularϕSampling` with ``N_\vartheta`` samples along ``\vartheta`` by either specifying ``J_\vartheta = 2 N_\vartheta`` or ``J_\vartheta = 2 N_\vartheta + 1``. The difference is that ``\Delta \vartheta`` is different for the two choices, such that the last ``\vartheta``-sample will either coincide with ``\vartheta_{N_\vartheta} = \pi`` or with ``\vartheta_{N_\vartheta} = \pi - \Delta \vartheta /2``.   
---

!!! note
    A `RegularθRegularϕSampling` has always multiple samples at the North Pole (one sample per ``\varphi``-value). However, only `RegularθRegularϕSampling`s with even values for ``J_\vartheta`` also feature samples at the South Pole. Odd values for ``J_\vartheta`` lead to no samples at the South Pole.
---

While local and even global interpolation is easy and efficient with a `RegularθRegularϕSampling`, the drawback is that evaluating spherical integrals is more costly as compared to a `GaussLegendreθRegularϕSampling`.

## [`GaussLegendreθRegularϕSampling`](@id gausslegendresampling)
Objects with a `GaussLegendreθRegularϕSampling` are sampled along ``\varphi`` with ``N_\varphi`` equally distributed samples  in ``0\leq \varphi < 2\pi`` and ``N_\vartheta`` ``\vartheta``-samples on a Gauß-Legendre based grid with ``0< \vartheta < \pi``, according to
- ``\vartheta_k=\text{arccos}(-x_k)``
where ``x_k`` are the roots of the ``N_\vartheta``th Legendre polynomial
- ``\varphi_n = (n-1)\, \Delta\varphi`` with ``n=1,\ldots,N_\varphi``

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

The two different `SphereSamplingStrategy`s are illustrated below for reference.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/compare_sampling_darkmode.svg" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/compare_sampling_darkmode.svg" width="750" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Illustration of different SphereSamplingStrategys. Grid points are marked by orange dots and the blue lines correspond to a RegularθRegularϕSampling for reference.  <br/>  <br/>  
    Left: RegularθRegularϕSampling with Jϕ = Jθ = 16, Right: GaussLegendreθRegularϕSampling with Nϕ = 16 and Nθ=8. 
  </figcaption>
</figure>
<br/>
```


[^1]: or multiple matrices, if the function ``\bm{f}(\vartheta, \varphi)`` is a vector function containing e.g., data for two polarizations

[^2]: I. Bogaert, Iteration-free computation of Gauss–Legendre quadrature nodes and weights, SIAM J.ournal on Scientific Computing, Vol. 36, No. 3, pp. A1008–A1026, 2014 [DOI. 10.1137/140954969](https://doi.org/10.1137/140954969)
