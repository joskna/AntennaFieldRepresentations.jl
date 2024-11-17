# Fast Computation of Field Interactions with the Multilevel Fast Multipole Method
In usual measurement scenarios, the probe antennas are much smaller than the antenna under test. The geometry of an exemplary measurement scenario is depicted below.
```@raw html
<figure>
<picture>
  <source  srcset="../assets/mlfmm_scenario.png" width="800">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Geometrical situation in a typical antenna measurement scenario.
  </figcaption>
</figure>
<br/>
```
The Multilevel Fast Multipole Method can be used to accelerate the calculation of the transmission between a large antenna under test and many probe antennas at different measurement positions. 
The algorithm for the computation of the antenna interactions is based on the [plane-wave representation of the ``S_{21}`` parameter between two antennas](@ref planewave_received)
```math
S_{21}=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\hat{\bm{F}}_1\left(\hat{\bm{k}} \right)\cdot
T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}} \cdot \hat{\bm{R}})\,
\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\, ,
```
where ``\hat{\bm{F}}_1\left(\hat{\bm{k}} \right)`` is the normalized far-field pattern of the transmitting antenna 1, ``\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)`` is the inverted far-field pattern of the receiving antenna 2, and
```math
T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}} \cdot \hat{\bm{R}})
=
-\sum \limits_{\ell' =0}^{L}\left(-\mathrm{j}\right)^{\ell'} \,
\left(2\ell'+1\right) \, 
\mathrm{h}_{\ell'}^{(2)}(k_0 \lvert \bm{R} \rvert)
\mathrm{P}_\ell'(\hat{\bm{k}} \cdot \hat{\bm{R}})
```
is the [plane-wave translation operator](@ref planewave_translation).

The algorithm for the fast computation of many ``S_{21}``-parameters in the described configuration has three stages.

In the aggregation stage, by using interpolation and phase shifts, the far-field pattern of the antenna under test is computed by combining the far-field patterns of many smaller subsets of the antenna.

The far-field pattern of the antenna under test is thereafter translated into a plane wave spectrum which can represent the incident fields for many probe positions simultaneously in a well-defined observation region.

Finally in the disaggregation step, the sampling rate of the ``k``-space representation of the far-field patterns of the antenna under test and the probe antennas is reduced in an anterpolation step before the integrals are evaluated at the low sampling density efficiently. Without the anterpolation step, all integrals would have to be evaluated with a huge sampling rate in the ``k``-space.

## [Aggregation of Far-Field Patterns](@id mlfmm_aggregation)
One of the key points for the performance of the algorithm is the fact that the far-field pattern of a large antenna can be composed from the individual far-field patterns of disjunct subsets of its equivalent current representations. The geometric situation of this so-called aggregation step is depicted in the figure below.
```@raw html
<figure>
<picture>
  <source  srcset="../assets/aggregate.png" width="500">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic representation of the aggregation step. The far-field pattern of the complete antenna under test (denoted by the large orange sphere) is composed of far-field patterns of subsets of the current densety (denoted by the smaller red spheres).
  </figcaption>
</figure>
<br/>
```

The subsets of the currents  each have a signifcantly smaller geometrical size than the current distribution for the complete antenna under test. Therefore, fewer plane-wave samples are sufficient to represent their individual far-field patterns. As a consequence, the computation of the (less densely sampled) far-field patterns from the subset currents is more efficient than the calculation of the (more densely sampled) far field pattern of the complete antenna under test.

The less densely sampled far-field patterns of the subset currents is interpolated to match the sampling density which is necessary to represent the complete far-field pattern of the antenna under test as indicated in the figure below.
```@raw html
<figure>
<picture>
  <source  srcset="../assets/interpolate.png" width="500">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic representation of the interpolation step.
  </figcaption>
</figure>
<br/>
```

The plane-wave samples of the densely sampled far-field patterns are then phase-shifted as ``\hat{\bm{F}}(\hat{\bm{k}}) \rightarrow \hat{\bm{F}}(\hat{\bm{k}}) \, \mathrm{e}^{\, \mathrm{j} \bm{k}\cdot \bm{R}}`` such that they are defined with respect to the same expansion center. 
```@raw html
<figure>
<picture>
  <source  srcset="../assets/phaseshift.png" width="500">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic representation of the phase-shifting step.
  </figcaption>
</figure>
<br/>
```

Finally, all interpolated and phase-shifted far-field patterns of the subset currents can be added together to form the far-field pattern of the complete antenna under test.



### Interpolation of Far-Field Patterns
Depending on the [sampling strategy of the far-field patterns](@ref planewave_sampling), the plane-wave samples of the densely sampled far-field pattern can be found from the plane-wave samples of the less densely sampled far-field pattern by subsequent 1D-interpolations along the ``\vartheta``- and ``\varphi``-component, respectively. 
Suitable techniques for 1D-interpolation are the trigonometric interpolation[^1] and piecewise polynomial interpolation[^2] [^3].

Trigonometric interpolation is particularly useful for equidistantly sampled periodic data, where the Fast Fourier Transform can accelerate the interpolation process. In this case, trigonometric interpolation is an exact interpolation technique (in exact arithmetic) with a computational complexity of ``\mathcal{O}(n \log n)`` to interpolate ``n`` points in a 1D periodic data set.  

If computational complexity is crucial, piecewise polynomial interpolation comes into play. A (low-order) polynomial is fitted into a small subset of sampling points around the point we wish to interpolate. The piecewise polynomial interpolation is not exact, opposed to the trigonometric interpolation, but the error can be controlled by the order of the polynomial fit and/or the distance between the sampling points. The computational complexity of the piecewise polynomial interpolation is ``\mathcal{O}(K\,n)`` to interpolate ``n`` points in a 1D data set, where ``K`` is the fixed order of the polynomial fit. Notice that the order of the polynomial does not depend on the overall length of the data set. Furthermore, the sampling is not required to be neither equidistantly spaced nor periodic (although the Runge phenomenon may obstruct accurate interpolation near the end points of the sampling for most sampling strategies [^3]). 

## Translation of the Antenna Under Test Far-Field Pattern into a Plane Wave Spectrum in Observation Region
After the far-field pattern of the antenna under test has been determined, we can calculate an incident plane wave spectrum for the observation region by applying the [translation operator](@ref planewave_translation). Ideally, the observation region is large and includes many probe antenna positions at once. In this case, the obtained incident plane wave spectrum is valid for calculating the received signals of all the probes within the observation region.

The valid region for the incident plane-wave representation is ultimately dictated by the contained number of spherical modes. With mire spherical modes in the field representation, the fields of the incident plane-wave spectrum can be correctly evaluated in a larger observation region. The drawback is that this will also neccesitate an increase  of the effective number of spherical modes contained in the antenna far-field pattern before the translation. Due to numerical instabilities of the higher-order spherical modes near the origin, we must also be further away from the antenna to evaluate the fields accurately in this case. 

As a rule of thumb, the center points of the antenna under test region and the observation region must be separated by at least ``2 \ell_{\mathrm{max}} / k_0 `` when ``\ell_{\mathrm{max}}`` is the maximum effective mode order for the ``\ell``-index of the spherical mode content in either of the field representations. For an accurate evaluation of the desired ``S_{21}-integral``
```math
S_{21}=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\hat{\bm{F}}_1\left(\hat{\bm{k}} \right)\cdot
T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}} \cdot \hat{\bm{R}})\,
\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\, ,
```
the summation in the translation operator must be carried out for all terms up to ``L=2 \ell_{\mathrm{max}}`` and the far-field patterns ``\hat{\bm{F}}_1\left(\hat{\bm{k}} \right)`` and  ``\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)`` should be sampled such that spherical mode orders up to ``\ell=2 \ell_{\mathrm{max}}`` can be resolved.

## Disaggregation of Plane-Wave Spectra
Consider the situation after the translation step as depicted in the figure below.
```@raw html
<figure>
<picture>
  <source  srcset="../assets/anterpolate.png" width="500">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Situation after the translation step.
  </figcaption>
</figure>
<br/>
```
The antenna under test field is represented by an incident plane-wave spectrum `` \hat{\bm{P}}_{\mathrm{inc}}\left(\hat{\bm{k}} \right) = \hat{\bm{F}}_1\left(\hat{\bm{k}} \right)T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}} \cdot \hat{\bm{R}}) `` in the observation region around the expansion center ``\bm{r}_0`` (denoted by the large green sphere). The probe antenna's far-field pattern is known with an expansion center at ``\bm{r}_m``. Furthermore, since the probe antenna is small, its far-field pattern has a low-order spherical mode content (thus the representation by a small blue sphere). To evaluate the interaction integral, we need to shift the probe pattern to the incident plane-wave expansion center ``\bm{r}_0`` and ensure that its ``k``-space representation has the same sampling as the incident plane-wave spectrum `` \hat{\bm{P}}_{\mathrm{inc}}\left(\hat{\bm{k}} \right) ``.

We have
```math
S_{21}=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\hat{\bm{P}}_{\mathrm{inc}}\left(\hat{\bm{k}} \right) 
\cdot
\mathrm{e}^{\,\mathrm{j}\bm{k} \cdot (\bm{r}_0-\bm{r}_m)}\,
\mathrm{interpolate}\left[
\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)
\right]
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\,, 
``` 
where we bring the probe far-field pattern to the correct sampling density by the interpolation ``\mathrm{interpolate}\left[\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)\right]`` before we phase-shift each far-field component by ``\mathrm{e}^{\,\mathrm{j}\bm{k} \cdot (\bm{r}_0-\bm{r}_m)}`` to correspond to the correct expansion center at ``\bm{r}_0``.

Let ``\text{anterpolate}[\cdot]`` be the [adjoint operation](https://en.wikipedia.org/wiki/Hermitian_adjoint) to the ``\text{interpolate}[\cdot]`` operation. Then we have
```math
S_{21}=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\mathrm{anterpolate}\left[
\hat{\bm{P}}_{\mathrm{inc}}\left(\hat{\bm{k}} \right)
\,
\mathrm{e}^{\,\mathrm{j}\bm{k} \cdot (\bm{r}_0-\bm{r}_m)}
\right] 
\cdot
\hat{\bm{F}}_2\left(- \hat{\bm{k}} \right)
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\,, 
``` 
i.e., the anterpolation operation applied to the (phase-shifted) incident plane-wave spectrum. If the interpolation algorithm can be represented by a matrix multiplication, the anterpolation algorithm is simply a multiplication by the Hermitian adjoint of the interpolation matrix. The result of the anterpolation operation is a plane-wave spectrum centered at the probe position ``\bm{r}_m`` matching the low sampling density of the original probe far-field pattern. In fact, the anterpolation step can be regarded as downsampling of the incident plane-wave spectrum with an inherent spherical-mode low-pass filter to avoid aliasing. Due to the neccesary spherical mode-filter, we cannot simply down-sample the incident plane-wave spectrum to match the probe sampling, but we must carefully implement the adjoint operator to our interpolation algorithm.

## [References](@id mlfmm_refs)
[^1]: Wikipedia: [Trigonometric interpolation](https://en.wikipedia.org/wiki/Trigonometric_interpolation)
[^2]: Wikipedia: [Polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation)
[^3]: J.-P. Berrut, L. N. Trefethen "Barycentric Lagrange Interpolation", SIAM Review, Vol. 46, No. 3, pp. 501-517, 2004 [DOI. 10.1137/S0036144502417715](https://doi.org/10.1137/S0036144502417715)