# Plane Wave Representations of Electromagnetic Fields


## Sampling of Plane Wave Representations of Electromagnetic Fields

Far field patterns ``\bm{F}(\vartheta, \varphi)`` and plane wave spectra ``\bm{P}(\vartheta, \varphi)`` are representated by a discrete set of samples in the ``\vartheta, \varphi``-domain. 

The samples should be chosen such that
- integration over the Ewald sphere is fast and accurate
- local interpolation is sufficently accurate (local interpolation is necessary for efficiency, e.g., with multipole based methods)
- the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` has small computational cost (ideally should not rely on interpolation)
- redundancy in the representation is kept minimal

The plane wave representation in this package have an equiangular sampling steps along ``\varphi`` with ``0\leq \varphi < 2\pi`` and distribute the ``\vartheta``-samples on a Gauß-Legendre based grid with ``0< \vartheta < \pi``, according to 
```math
\vartheta_k=\text{arccos}(-x_k)
```
where ``x_k`` are the roots of the ``n``th Legendre polynomial. To represent any far-field pattern or plane-ave spectrum up to a spherical mode order of ``L``, the number of ``\vartheta``- and ``\varphi``-samplesis 
```math
N_\vartheta=L+1
```
```math
N_\varphi=2L+2
```
Although the ``\varphi``-sampling is slightly redundant (``2L+1`` samples would suffice) the even number of  ``\varphi``-samples ensures that for every sampling point there is another sampling point exactly in the oppsoite direction. This has obvious benefits for the conversion  ``\bm{F}(\hat{\bm{k}}) \rightarrow \bm{F}(-\hat{\bm{k}})`` as well as for pattern interpolation, when the interpolation scheme requires samples at either side of the coordinate poles at ``\vartheta\in \{0,\pi\}``.
The figures below illustrate the regular sampling points and the Lagrange sampling points utilized by `PlaneWaveRepresentation`s. 


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

## Integration
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

## Interpolation

Polynomial interpolation finds a polynom ``p_n(x)`` of order ``n`` such that 
```math
p_n(x_i)=f(x_i) \qquad i=1,\ldots,n \, .
```
The polynomial can b expressed as[^2]

```math
p_n(x)=\sum \limits_{i=1}^n f(x_i) w_i(x)
```
with 
```math
w_i(x) =\dfrac{\prod \limits_{j=1}^n (x-x_j)}{\prod \limits_{j=1, j\neq i}^n (x_i-x_j)}\, .
```
Fairly accurate and efficient local interpolation can be achieved with the barycentric Lagrange interpolation of moderate order applied to periodic extentions of the pattern functions along circles with fixed ``\varphi`` and ``\vartheta`` (in particular the Legendre points are known to lead to well-conditioned polynomial interpolations [^3], relatively immune against the Runge phenomenon which is commonly observed with equispaced sampling points on limited sampling domains). 


## [References](@id refs)
[^1]: I. Bogaert, Iteration-free computation of Gauss–Legendre quadrature nodes and weights, SIAM J.ournal on Scientific Computing, Vol. 36, No. 3, pp. A1008–A1026, 2014 [DOI. 10.1137/140954969](https://doi.org/10.1137/140954969)
[^2]: Wikipedia: [Polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation)
[^3]: J.-P. Berrut, L. N. Trefethen "Barycentric Lagrange Interpolation", SIAM Review, Vol. 46, No. 3, pp. 501-517, 2004 [DOI. 10.1137/S0036144502417715](https://doi.org/10.1137/S0036144502417715) 
