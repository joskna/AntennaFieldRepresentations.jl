# [Interpolation of Spherically Sampled Data] (@id interpolation_theory)
Consider a function ``f(\vartheta, \varphi)`` defined on the surface of a sphere, which is sampled at a discrete set of sampling points ``(\vartheta_m, \varphi_n)`` with ``0\leq \vartheta_m \leq \pi``, ``0 \leq \varphi_n < 2\pi``, ``m=1,\ldots, N_\vartheta``, ``n=1,\ldots, N_\varphi``. The sampling points are illustrated below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/sampling_regular_interp_darkmode.svg" width="200">
  <source media="(prefers-color-scheme: light)" srcset="../assets/sampling_regular_interp_darkmode.svg" width="200" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Illustration of the sampling points for the sampled function on the sphere. 
  </figcaption>
</figure>
<br/>
```

In many circumstances, we need to evaluate the function ``f(\vartheta, \varphi)`` at angles ``(\vartheta, \varphi)`` which do not coincide with any of the stored data points ``f(\vartheta_m, \varphi_n)``. In this case, we must perform an interpolation of the data to the desired location. 
As long as the data is sampled on closed rings with constant ``\vartheta`` or constant ``\varphi`` [^1] , respectively, one can perform the interpolation as a sequence of 1D interpolations along ``\vartheta`` and ``\varphi``. Two closed rings with constant ``\vartheta`` or constant ``\varphi`` are illustrated below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/sampling_regular_interp_darkmod_highlighted.svg" width="200">
  <source media="(prefers-color-scheme: light)" srcset="../assets/sampling_regular_interp_darkmod_highlighted.svg" width="200" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Illustration of sampled data on closed rings with constant Theta (green) or constant Phi (red). 
  </figcaption>
</figure>
<br/>
```
In the following, we will only outline the 1D interpolation procedures for ``2\pi``-periodic functions.
We will distinguish between *global interpolation*, i.e., an interpolation scheme which uses all available data points along a closed ring (i.e., a ring with either a fixed ``\vartheta``-coordinate or a fixed ``\varphi``-coordinate), and *local interpolation*, i.e. an interpolation scheme which utlizes only a few sampling points in the neighbourhood of the desired function input.

## Global Interpolation
Global interpolation leverages on the fact that any band limited ``2\pi``-periodic function ``f(\vartheta)`` may be expressed via

```math
f(\vartheta)
=
\sum
\limits_{k=-B}^{B}
c_k \, \mathrm{e}^{\,\mathrm{j} k \theta}\, .
```

If we have at least ``J_\vartheta = 2B+1`` sampling points on the ring, we can, in principle, reconstruct all coefficients ``c_k``. From the reconstructed coefficients ``c_k`` one could calculate the function ``f(\vartheta)`` at any arbitrary angle ``\vartheta`` using the equation above.

If the ``J_\vartheta`` input sampling points are equidistantly spaced, we can use an FFT to calculate the coefficients ``c_k`` with ``\mathcal{O}(J_\vartheta \log J_\vartheta)`` computational complexity. The calculation of the new function value from the coefficients ``c_k`` requires an additional ``\mathcal{O}(J_\vartheta)`` operations, such that the overall interpolation process is completely governed by the computional cost of the FFT.

!!! note
    The computational complexity of the interpolation of **equidistantly sampled periodic data** is ``\mathcal{O}(J_\vartheta \log J_\vartheta)`` because it is determined by the cost of an FFT. 
---

If we further want to interpolate the function ``f(\vartheta)`` at ``J_\vartheta' `` equidistantly spaced new points, we can fill zeros into the array containing the coefficients ``c_k`` such thath the length of the array matches the desired number of sampling points.
All of the new function values can be obtained at once from the zero-padded vector of coefficients ``c_k`` with a single IFFT, such that the resulting sampling points  can be computed in ``\mathcal{O}(J_\vartheta' \log J_\vartheta')`` computational complexity.

!!! note
    If we want to interpolate the function ``f(\vartheta)`` at ``J_\vartheta' `` **equidistantly spaced new points** from sampled data at 
    ``J_\vartheta `` **equidistantly spaced old sample** points ``f(\vartheta_n)``, the computational cost is ``\mathcal{O}(J_\vartheta \log J_\vartheta) + \mathcal{O}(J_\vartheta' \log J_\vartheta')`` 
---

However, depending on the `SphereSamplingStrategy`, the sampling points may or may not be equidistantly spaced. For Non-equidistantly spaced sampling grids, we can use the Barycentric Trigonometric Interpolation formula[^2]
```math
f(\vartheta) = 
\dfrac{\sum \limits_{k=1}^{J_\vartheta} \dfrac{w_k}{\sin\left(\frac{1}{2} (\vartheta -\vartheta_k) \right)}\, f(\vartheta_k)}{\sum \limits_{k=1}^{J_\vartheta} \dfrac{w_k}{\sin\left(\frac{1}{2} (\vartheta -\vartheta_k )\right)}}\, ,
```
where 
```math
w_k= 
\dfrac{1}{\prod \limits_{i\neq k} \sin\left(\frac{1}{2} (\vartheta_k -\vartheta_i )\right)};
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~i =1, \ldots, N_\vartheta\,.
```
This interpolation formula is also valid for non-equidistantly spaced sampling points ``\vartheta_k``. 

!!! warning
    This interpolation formula is only valid for odd ``J_\vartheta``. If we have an even number of irregularly spaced sampling points ``J_\vartheta``, it is recommended to discard one of the sampling points.
---

!!! note
    Note that the barycentric trigonometrc interpolation formula takes the form ``f(\vartheta) = \sum_{k=1}^{J_\vartheta} \alpha_k\, f(\vartheta_k)``. With predefined interpolation weights ``\alpha_k`` the computational cost of evaluating the global interpolation is ``\mathcal{O}(J_\vartheta)``, i.e., the computational cost depends on the size of the input data.
---

!!! note
    If we want to interpolate the function ``f(\vartheta)`` at ``J_\vartheta' `` new points, we need to evaluate the barycentric trigonometrc interpolation formula ``J_\vartheta' `` times, leading to a computational cost of ``\mathcal{O}(J_\vartheta \,J_\vartheta')``. 
---

For bandlimited functions ``f(\vartheta)``, global interpolation should be exact (neglecting rounding errors).

## Local Interpolation
In certain situations, it is not feasible to have the computational cost of the interpolation to be dependent on the size of the input data. In these cases, we rely on local (polynomial) interpolation [^3] [^4]. 

While global polynomial interpolation fits a polynomial to the entire data on the ring (thus, taking into account all available data points), local interpolation fits a specified order polynomial (fourth order, fifth order, sixth order, ... etc.) using only points in the adjacient neighbourhood of the desired function input ``\vartheta_{\mathrm{new}}``. For example, a sixth-order local interpolation would fit a sixth order polynomial through the three sample points before and after ``\vartheta_{\mathrm{new}}``. In general, a ``K``-order local polynomial will utilize the adjacient ``K`` data points to interpolate the function ``f(\vartheta)`` at the new input ``\vartheta``. The situation is illustrated below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/local_interp_1Dcircle.svg" width="300">
  <source media="(prefers-color-scheme: light)" srcset="../assets/local_interp_1Dcircle.svg" width="300" >
  <img alt="" src="" width="300">
</picture>

  <figcaption>
    To interpolate the function at the green location, from all available data points (faint orange) only the highlighted sample points (solid red) are utilized for the interpolating polynomial.
  </figcaption>
</figure>
<br/
```

We have [^4]
```math
f(\vartheta) \approx
\dfrac{\sum \limits_{k=1}^{K} \dfrac{w_k}{\vartheta - \vartheta_k}\, f(\vartheta_k)}{\sum \limits_{k=1}^{K} \dfrac{w_k}{\vartheta - \theta_k}}
```
where ``\vartheta_k`` are the ``K`` old sampling points in the neighbourhood of ``\vartheta`` and
```math
w_k= 
\dfrac{1}{\prod \limits_{i\neq k}  (\vartheta_k -\vartheta_i )};
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~i =1, \ldots, K\,.
```

This interpolation formula is also valid for equidistantly and non-equidistantly spaced sampling points ``\vartheta_k``.

!!! note
    Note that the local interpolation formula takes the form ``f(\vartheta) = \sum_{k=1}^{K} \alpha_k\, f(\vartheta_k)``. With predefined interpolation weights ``\alpha_k`` the computational cost of evaluating the global interpolation is ``\mathcal{O}(K)``, i.e., the computational cost *does not dependd on the size of the input data*.
---

!!! note
    If we want to interpolate the function ``f(\vartheta)`` at ``J_\vartheta' `` new points, we need to evaluate the local interpolation formula ``J_\vartheta' `` times, leading to a computational cost of ``\mathcal{O}(K \,J_\vartheta')``. 
---

Although local interpolation is not exact for bandlimited functions ``f(\vartheta)``, the approximation error can be controlled by 
- oversampling of the original data (i.e., the data is originally stored with more samplng points than required by the Nyquist critereon for recovery of the underlying bandlimited function)
- increasing the degree of the interpolating polynomial

[^1]: For an equidistant sampling along the ``\varphi``-direction, this requires an even number of ``\varphi``-samples ``N_\varphi`` such that the ``\vartheta``-rings can be interpolated over the poles.

[^2]: Berrut, Jean -Paul. "Baryzentrische Formeln zur trigonometrischen Interpolation (I)." Zeitschrift für angewandte Mathematik und Physik ZAMP 35 (1984): 91-105 [DOI. 10.1007/BF00945179](https://doi.org/10.1007/BF00945179)

[^3]: Trefethen, L. "Six myths of polynomial interpolation and quadrature." (2011). [University of Oxford](https://people.maths.ox.ac.uk/trefethen/mythspaper.pdf)

[^4]: Berrut, Jean-Paul, and Lloyd N. Trefethen. "Barycentric lagrange interpolation." SIAM review 46.3 (2004): 501-517. [DOI. 10.1137/S003614450241771](https://epubs.siam.org/doi/abs/10.1137/S0036144502417715)
