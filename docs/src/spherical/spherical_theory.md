# Spherical Vector-Wave Expansion of General Time Harmonic Electromangentic Fields
Away from sources, the time harmonic electric and magnetic fields must fulfill the homogeneous curl-curl-equation 
```math
\nabla \times \nabla \times \bm{E} - k_0^2 \bm{E} =\bm{0}\, .
```

Arbitrary source-free solutions of Maxwell's equations can be expanded into a series of spherical vector-wave functions as (see pp. 9ff of [^1])
```math
\bm{E}(r, \vartheta, \varphi) =
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(c)}\,
\bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi)
```
```math
\bm{H}(r, \vartheta, \varphi) =\mathrm{j}\,
\dfrac{k_0}{\sqrt{Z_{\mathrm{F}}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(c)}\,
\bm{F}_{3-s,\ell m}^{(c)}(r,\vartheta, \varphi)
```  
where ``Z_{\mathrm{F}}\sqrt{\mu_0/ \varepsilon_{0}}\approx 376.730\,313\,669\, \mathrm{\Omega}`` is the impedance of free space and ``\bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi)`` are the spherical vector wave mode functions, which are solutions to the homogeneous curl-curl equation (except in the origin of the coordinate system, where some types of vector mode functions become singular). Due to only two of the radial dependencies being linearly independent, it is sufficient to let the sum over ``c`` to be an arbitrary pair of two iindices from the set ``c \in \{1,2,3,4\}``. Usually one will choose the pair ``c=1`` and ``c=4``, where only the ``c=4``-type modes are needed to describe purely radiated fields and only the ``c=1``-type modes are needed to describe purely incident fields.

In the numerical implementation, the coefficients ``\alpha_{s \ell m}`` are stored in a vector, where the triple index
``s \ell m`` is mapped to a single consecutively running index ``j`` by the function `sℓm_to_j(s,ℓ,m)`. The inverse map from a single index to the triple index ``s \ell m`` is implemented in the function `j_to_sℓm(j)`.

The vector of coefficients completely defines the fields of a corresponding radiating or incident spherical wave expansion. 



## [Definition of Spherical Vector-Mode Functions](@id spherical_definition)

Explicit expressions for the spherical vector-mode functions are (see p. 13 of [^2])
```math
\bm{F}_{1\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}
\left(
\dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \, \dfrac{\mathrm{j}m \,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\sin \vartheta} \bm{e}_\vartheta
-\dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \, \dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta} \bm{e}_\varphi
\right)
```
```math
\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}
\left(
\dfrac{\sqrt{\ell(\ell+1)}}{k_0r}\, z_\ell^{(c)}(k_0r)\, \overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) \bm{e}_r
\right.
\phantom{\dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta}~~~~~}
\\
\phantom{\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}}
+ \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \right\}\, \dfrac{\mathrm{d}\,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\mathrm{d} \vartheta} \bm{e}_\vartheta
\phantom{~~~~~}
\\ \left.
\phantom{\bm{F}_{2\,\ell m}^{(c)}(r,\vartheta, \varphi)
=
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2}\pi}}
+ \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \dfrac{z_\ell^{(c)}(k_0r)}{\sqrt{\ell(\ell+1)}} \right\}\, \dfrac{\mathrm{j}m \,\overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) }{\sin \vartheta} \bm{e}_\varphi\
\right)\, .
```


To ensure that the vector-mode functions ``\bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi)`` are solutions to the homogeneous curl-curl equation, they are constructed from solutions of the homogeneous scalar Helmholtz equation as (see, e.g., page 393 of [^3])
```math
\bm{M}_{\ell m}^{(c)}= \nabla f_{\ell m}^{(c)} \times \bm{e}_r \overset{!}{=} \bm{F}_{1\,\ell m}^{(c)}
```
and
```math
\bm{N}_{\ell m}^{(c)}= \dfrac{1}{k_0} \nabla \times \bm{M}_{\ell m}^{(c)} \overset{!}{=} \bm{F}_{2\,\ell m}^{(c)}
``` 
where ``\bm{e}_r`` is the unit vector in radial direction.


General solution to the scalar Helmholtz equation 
```math
\Delta f +k_0^2 f=0\, .
```
can be constructed as superpositions of functions from the set
```math
f_{s\ell m}^{(c)}(r, \vartheta, \varphi)
=
\begin{cases}
\dfrac{\mathrm{j}_{\ell}(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=1\\[1.5em]
\dfrac{\mathrm{\mathrm{n}}_\ell(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=2\\[1.5em]
\dfrac{\mathrm{\mathrm{h}}_\ell^{(1)}(k_0 r)}{\sqrt{\ell(\ell+1)}}
\, \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=3\\[1.5em]
\dfrac{\mathrm{\mathrm{h}}_\ell^{(2)}(k_0 r)}{\sqrt{\ell(\ell+1)}}\,
 \mathrm{Y}_{\ell m}(\vartheta, \varphi) & c=4
\end{cases} 
```
where 
- ``\mathrm{j}_{\ell}(k_0 r)`` are spherical Bessel functions which represent the radial behavior standing waves (or *incident* waves)
- ``\mathrm{n}_{\ell}(k_0 r)`` are spherical Neumann functions which usually do not represent a physically meaningful solution
- ``\mathrm{h}_\ell^{(1)}(k_0 r) = \mathrm{j}_{\ell}(k_0 r) +\mathrm{j}\, \mathrm{n}_{\ell}(k_0 r)`` are spherical Hankel functions of first kind which represent the radial behavior of waves coming from infinity being absorbed at the coordinate origin
- ``\mathrm{h}_\ell^{(2)}(k_0 r) = \mathrm{j}_{\ell}(k_0 r) -\mathrm{j}\, \mathrm{n}_{\ell}(k_0 r)`` are spherical Hankel functions of second kind which represent the radial behavior of waves traveling from the coordinate origin toward infinity
- ``\mathrm{Y}_{\ell m}(\vartheta, \varphi)`` are the spherical Harmonics
and ``k_0=2\pi/\lambda`` is the wavenumber defined by the frequency of the electromagnetic wave.



The spherical Harmonics are defined as (see §14.30 of [^1])
```math
\mathrm{Y}_{\ell m}(\vartheta, \varphi) 
=
\begin{cases}
\dfrac{1}{\sqrt{2\pi}}\, \overline{\mathrm{P}}_{\ell}^{m}(\cos \vartheta) \mathrm{e}^{\, \mathrm{j}m \varphi} & \text{for }|m|\leq \ell\\
0 & \text{esle}
\end{cases}
``` 
where the normalized associated Legendre functions ``\overline{\mathrm{P}}_{\ell}^{m}(x)`` are defined in ordinary Legendre polynomials ``\mathrm{P}_{\ell}(x)`` 
as (compare §14.6 and §14.7 of [^1])
```math
\overline{\mathrm{P}}_{\ell}^{m}(x)
=
\begin{cases}
\sqrt{\dfrac{(2\ell+1)(\ell-m!)}{2(\ell+m)!}}
(-1)^m
\,
(1-x^2)^{m/2}
\, \dfrac{\mathrm{d}^m \mathrm{P}_{\ell}(x)}{\mathrm{d}^m} & m\geq 0\\[2em]
(-1)^m  \overline{\mathrm{P}}_{\ell}^{-m}(x) & m<0\, .
\end{cases}
```

!!! note
    Different sign conventions are used by different authors when the normalized associated Legendre functions are defined. Conventions may differ by a sign-factor of ``(-1)^m`` (this factor is often called *Condon-Shortley phase factor*). We follow the sign convention of [^1], which includes the Condon-Shortley phase factor of ``(-1)^m`` in the definition of the normalized associated Legendre functions for positive ``m``. The main reference [^2] for our implementation does not include the phase factor in the definition for the associated Legendre functions, but includes an explicit ``(-1)^m``-term to the definition of the mode functions (which we don't do), such that in the end, our definition of the final mode functions coincides with the definition used in [^2].
---
!!! warning
    A straightforward implementation of the above definition for the normalized associated Legendre functions will fail already for moderately large mode orders ``\ell`` due to the involved factorials. A more sensible definition is based on recurrence relations[^4]
---

## [Orthogonality Relations of Spherical Vector-Mode Functions](@id spherical_orthogonality)
The tangential components of the spherical vector-mode functions are orthogonal on concentric spherical surfaces (centered at the coordinate origin). Concretely, they fulfill the orthogonality properties
```math
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
  \bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi) \bigg\rvert_{\mathrm{tan}}
\,\cdot\,
\bm{F}_{\sigma \lambda \mu}^{(\gamma)}(r,\vartheta, \varphi) \bigg\rvert_{\mathrm{tan}}
\, 
\sin \vartheta\,
\mathrm{d}\vartheta
\mathrm{d}\varphi
= 
\delta_{s\sigma}\, 
\delta_{\ell \lambda}\, 
\delta_{m, -\mu}\,
R_{s\ell}^{(c)}(k_0r)\,
R_{\sigma\\lambda}^{(\gamma)}(k_0r)\,,
``` 
where
```math
\delta_{\ell \lambda}
=
\begin{cases}
1 & \text{if } \ell= \lambda \\
0 & \text{if } \ell\neq \lambda
\end{cases}
```
is the Kronecker-Delta and
```math
R_{s\ell}^{(c)}(k_0r)
= \begin{cases}
z_{\ell}^{(c)}(k_0r) & \text{if } s=1 \\[1em]
\dfrac{1}{k_0r}\, \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{z_{\ell}^{(c)}(k_0r)\right\} & \text{if } s=2
\end{cases}
```
is a short-hand notation for the radial dependency of the corresponding spherical vector-mode function ``\bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi)``.

The orthogonality relations are very useful to find the coefficients ``\alpha_{s\ell m}^{(c)}`` in the spherical vector-wave expansion of a certain field which is known on a spherical surface.

## [Far-Field Expressions](@id spherical_farfield)
Due to the asymptotic behavior of the radial dependencies 
```math
\lim \limits_{k_0r \rightarrow \infty} \mathrm{h}_\ell^{(2)}(k_0r) = \mathrm{j}^{\ell+1} \dfrac{\mathrm{e}^{-\mathrm{j} k_0r}}{k_0r} 
```
```math
\lim \limits_{k_0r \rightarrow \infty}  \dfrac{1}{k_0r} \dfrac{\mathrm{d}}{\mathrm{d} k_0r} \left\{k_0 r \mathrm{h}_\ell^{(2)}(k_0r) \right\} = \mathrm{j}^{\ell} \dfrac{\mathrm{e}^{-\mathrm{j} k_0r}}{k_0r} 
```
for distances approaching infinity, it is convenient to define far-field expressions of the mode functions as
```math
\bm{K}_{1\,\ell m}^{(4)}(\vartheta, \varphi) =
\lim \limits_{k_0r \rightarrow \infty}
\dfrac{k_0r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{F}_{1\,\ell m}^{(4)}(r,\vartheta, \varphi)
```
```math
\bm{K}_{2\,\ell m}^{(4)}(\vartheta, \varphi) =
\lim \limits_{k_0r \rightarrow \infty}
\dfrac{k_0r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{F}_{2\,\ell m}^{(4)}(r,\vartheta, \varphi)\, .
```

Explicit expressions for the far-field mode functions are given as
```math
\bm{K}_{1\,\ell m}^{(4)}(\vartheta, \varphi) =
\dfrac{\mathrm{j}^{\ell+1}}{\sqrt{\ell(\ell+1)}}
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2\pi}}
\left(
\dfrac{\mathrm{j}m\, \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\sin \vartheta} \bm{e}_\vartheta
-
\dfrac{\mathrm{d} \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\mathrm{d} \vartheta} \bm{e}_\varphi  
\right)
```
```math
\bm{K}_{2\,\ell m}^{(4)}(\vartheta, \varphi) =
\dfrac{\mathrm{j}^{\ell}}{\sqrt{\ell(\ell+1)}}
\dfrac{\mathrm{e}^{\, \mathrm{j}m\varphi}}{\sqrt{2\pi}}
\left(
\dfrac{\mathrm{d} \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\mathrm{d} \vartheta} \bm{e}_\vartheta
+
\dfrac{\mathrm{j}m\, \overline{\mathrm{P}}_{\ell}^m(\cos \vartheta )}{\sin \vartheta}
 \bm{e}_\varphi  
\right)\, .
```
!!! note
    This definition of the far-field mode functions lacks the arbitrary factor of ``\sqrt{4\pi}`` which is introduced in [^2] .
---


Naturally, the far field ``\bm{F}(\vartheta, \varphi) = \lim \limits_{k_0r \rightarrow \infty} \dfrac{r}{\mathrm{e}^{-\mathrm{j} k_0r}}
\bm{E}(r,\vartheta, \varphi) `` can be calculated from a given set of expansion coefficients ``\alpha_{s\ell m}^{(4)}`` as
```math
\bm{F}(\vartheta, \varphi) =
\sqrt{Z_{\mathrm{F}}}
\sum \limits_{c}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(4)}\,
\bm{K_{s\ell m}^{(c)}}(\vartheta, \varphi)\, .
```

## Spherical Expansion Coefficients from a Current Density
The radiating spherical vector-wave expansion for the fields radiated by arbitrary electric and magnetic volume current densities ``\bm{J}(\bm{r})`` and ``\bm{M}(\bm{r})`` can be calculated by 
```math
\alpha_{s\ell m}^{(4)}=k_0\, (-1)^{m+1} 
\iiint 
\left[
\sqrt{Z_{\mathrm{F}}}\, \bm{F}_{s,\ell,-m}^{(1)} \cdot  \bm{J}(\bm{r})
-
\dfrac{\mathrm{j}}{\sqrt{Z_{\mathrm{F}}}}(\bm{r})\, \bm{F}_{3-s,\ell,-m}^{(1)}(\bm{r}) \cdot  \bm{M}(\bm{r})
\right] 
\mathrm{d}v \, .
```

The spherical vector-wave functions ``\bm{F}_{s\ell m}^{(1)}(\bm{r})`` of incident type which are involved in the above calculation have vanishingly small magnitudes for locations ``\bm{r}`` with distance to the coordinate origin ``r < \ell /k_0``. Thus, for current densities which are completely enclosed by a minimum sphere of radius ``r_{\mathrm{minsphere}}``, one may neglect all coefficients ``\alpha_{s\ell m}^{(4)}`` with ``\ell < k_0 r_{\mathrm{minsphere}}`` (plus a small buffer depending on the desired accuracy), as long as one is not interested in evaluating the expanded electromagnetic fields extremely close to the original currents. 

!!! note
    The physical size of the current distribution dictates the maximum degree ``\ell`` of the spherical wave expansion to accurately characterize the radiated fields.
---

## [Spherical Expansion Coefficients by Using Orthogonality](@id spherical_coefficients_orthogonality)
Often, we are interested in finding the spherical vector-wave expansion of a purely incident or a purely radiated field. To be precise, one wants to find the expansion coefficients ``\alpha_{s\ell m}^{(c)}`` (with ``c=1`` for a purely incident field and ``c=4`` for a purely radiated field) to fulfill the equality
```math
\bm{E}(r, \vartheta, \varphi) =
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(c)}\,
\bm{F}_{s\ell m}^{(c)}(r,\vartheta, \varphi)
```
for a given field ``\bm{E}(r, \vartheta, \varphi)``.

As indicated [previously](@ref spherical_orthogonality), one may leverage on the orthogonality properties of the spherical vector-wave functions to find the desired expansion coefficients ``\alpha_{s\ell m}^{(c)}``. We have
```math
\alpha_{s\ell m}^{(c)}
   =
   \dfrac{(-1)^{m}}{k_0\sqrt{Z_{\mathrm{F}}}\,R_{s\ell}^{(c)}(k_0r)\,
     R_{s\ell}^{(\gamma)}(k_0r)}
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
  \bm{E}(r,\vartheta, \varphi) \bigg\rvert_{\mathrm{tan}}
\,\cdot\,
\bm{F}_{s \ell ,-m}^{(\gamma)}(r,\vartheta, \varphi) \bigg\rvert_{\mathrm{tan}}
\, 
\sin \vartheta\,
\mathrm{d}\vartheta
\mathrm{d}\varphi\, .
```
The above formula for finding the spherical coefficients ``\alpha_{s\ell m}^{(c)}`` is valid at any radial distance ``r`` from the origin. In particular, the formula can also be used to find the coefficients ``\alpha_{s\ell m}^{(c)}`` from the radiated fields at far-field distance.

## [Number of Relevant Modes of an Incident Electromagnetic Field](@id spherical_incident_expansion)
In a source-free region around the coordinate origin, incident fields can be completely described in terms of vector spherical wave functions ``\bm{F}_{s\ell m}^{(1)}(\bm{r})`` of incident type by (accordingly for the H-field)
```math
\bm{E}(r, \vartheta, \varphi) =
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{\infty}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(1)}\,
\bm{F}_{s\ell m}^{(1)}(r,\vartheta, \varphi)\,.
```

The spherical vector-wave functions ``\bm{F}_{s\ell m}^{(1)}(\bm{r})`` of incident type which are involved in the above calculation have vanishingly small magnitudes for locations ``\bm{r}`` with distance to the coordinate origin ``r < \ell /k_0``. Thus, as long as one is interested only in the fields inside a spherical region around the origin with radius ``r_{\mathrm{observation}}``, it is sufficient to consider only modes with ``\ell < L = k_0 r_{\mathrm{observation}}`` (plus a small buffer depending on the desired accuracy). 
We have
```math
\bm{E}(r, \vartheta, \varphi)\,\bigg\lvert_{r< \ell/k_0} \approx
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{L}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(1)}\,
\bm{F}_{s\ell m}^{(1)}(r,\vartheta, \varphi) \,.
```

!!! note
    The physical size of the observation region around the coordinate origin dictates the maximum degree ``\ell`` of the spherical wave expansion to accurately characterize the incident fields.
---

## [Translation of Radiated Modes into Incident Modes in Different Coordinate System](@id spherical_translation)
The radiated electromagnetic fields generated by currents around the coordinate origin may be represented by a superposition of spherical vector-wave functions ``\bm{F}_{s\ell m}^{(4)}(r, \vartheta, \varphi)`` of radiated type. However, in a translated coordinate system (denoted by the primed coordinates ``r', \vartheta', \varphi' ``), the sources are far away from the new coordinate origin and the electromagnetic fields around the new origin may be represented in terms of spherical vector-wave functions ``\bm{F}_{s\ell m}^{(1)}(r', \vartheta', \varphi')`` of incident type with respect to the new coordinate origin. The geometrical situation is depicted in the figure below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/SphericalTranslation_dark.svg" width="500">
  <source media="(prefers-color-scheme: light)" srcset="../assets/SphericalTranslation.svg" width="500" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic of the geometry for the translation of radiated modes into incident modes in a different coordinate system.
  </figcaption>
</figure>
<br/>
```

The expansion coefficients ``\alpha_{\sigma \lambda \mu}^{(1)}`` of the incident field in the new coordinate system can be calculated from the expansion coefficients ``\alpha_{s\ell m}^{(4)}`` of the radiated field in the original coordinate system by
```math
   \alpha_{\sigma \lambda \mu}^{(1)}
   = 
   \sum \limits_{s,\ell, m}
   \mathcal{T}_{s\ell m}^{\sigma \lambda \mu}(\bm{R})\,
    \alpha_{s\ell m}^{(4)}
```
where the translation coefficient ``\mathcal{T}_{s\ell m}^{\sigma \lambda \mu}(\bm{R})`` is given by
```math
\mathcal{T}_{s\ell m}^{\sigma \lambda \mu}(\bm{R})
=
-
  \sum \limits_{\ell'=0}^\infty
  (-1)^{\sigma + \lambda +\mu}
  \int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
  (-\mathrm{j})^{\ell'}
  \left(
  2\ell'+1
  \right)
    {\mathrm{h}_{\ell'}^{(2)}}\left({k_0\lvert{\bm{R}}\rvert}\right)\,
  {\mathrm{P}_{\ell'}}\left({\hat{\bm{k}}\left(\vartheta, \varphi\right) \cdot\hat{\bm{R}}}\right)\,  
  {\bm{K}_{s\ell m}^{(4)}}\left(\vartheta, \varphi\right)
  \cdot
  {\bm{K}_{\sigma, \lambda,- \mu}^{(4)}}\left(\vartheta, \varphi\right) 
  \, \sin \vartheta 
  \mathrm{d}\vartheta\, \mathrm{d}\varphi\, ,
```
where ``\hat{\bm{k}}\left(\vartheta, \varphi\right)`` denotes the radial unit vector defined by the angles ``\vartheta, \varphi`` and ``\hat{\bm{R}}`` denotes the unit vector in the direction of the vector from the original coordinate origin to the new coordinate origin.

## Received Signal
Let ``\hat{\alpha}_{s\ell m}^{(c)}`` be the normalized expansion coefficients for the fields radiated by an antenna such that ``\alpha_{s\ell m}^{(c)}=  a\,\hat{\alpha}_{s\ell m}^{(c)}`` are the actual expansion coefficients when the antenna's transmit port is excited by a signal with wave amplitude ``a \in \mathbb{C} \, \sqrt{\mathrm{W}}``.
To determine the received signal ``b \in \mathbb{C} \, \sqrt{\mathrm{W}}`` of this antenna in receive mode under a certain illumination of an incident field, it is convenient to represent the incident field via its spherical vector wave expansion
```math
\bm{E}(r, \vartheta, \varphi)=
k_0 \, \sqrt{Z_{\mathrm{F}}}
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{L}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(1)}\,
\bm{F}_{s\ell m}^{(1)}(r,\vartheta, \varphi) \,.
```
In this case, the received signal is given by
```math
b=
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{L}
\sum \limits_{m=-\ell}^{\ell}
\alpha_{s\ell m}^{(1)}\,
\hat{\beta}_{s\ell m}, ,
```
where the antenna receive coefficients ``\hat{\beta}_{s\ell m}`` are found from the corresponding normalized antenna's transmit coefficients ``\hat{\alpha}_{s\ell m}^{(c)}`` via
```math
\hat{\beta}_{s\ell m}
=
\dfrac{(-1)^{m}}{2}
\hat{\alpha}_{s\ell,-m}^{(c)}\, .
```

Combining this result with the findings from the [section above](@ref spherical_translation), we can express the ``S_{21}`` parameter measured for the transmission between two antennas as
```math
S_{21}
=
\sum \limits_{\sigma=1}^2
\sum \limits_{\lambda=1}^{L}
\sum \limits_{\mu=-\lambda}^{\lambda}
\hat{\beta}_{\sigma \lambda \mu}^{\text{antenna2}}
\,
\left(
\sum \limits_{s=1}^2
\sum \limits_{\ell=1}^{L}
\sum \limits_{m=-\ell}^{\ell}
   \mathcal{T}_{s\ell m}^{\sigma \lambda \mu}(\bm{R})\,
    \hat{\alpha}_{s\ell m}^{(4),\text{antenna1}}
\right)
```
where ``\hat{\alpha}_{s\ell m}^{(4),\text{antenna1}}`` are the normalized transmit coefficients of antenna 1 (the transmit antenna), ``\hat{\beta}_{\sigma \lambda \mu}^{\text{antenna2}}`` are the receive coefficients of antenna 2 (the receive antenna) and ``\mathcal{T}_{s\ell m}^{\sigma \lambda \mu}(\bm{R})`` is the translation operator known from [previous sections](@ref spherical_translation) with ``\bm{R}`` denoting the vector separating the centers of the two antennas.
The geometric situation for the calculation of the ``S_{21}``-parameter is depicted below. 

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/S21Translation_dark.svg" width="500">
  <source media="(prefers-color-scheme: light)" srcset="../assets/S21Translation.svg" width="500" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The geometric situation for the calculation of the transmission between two antennas.
  </figcaption>
</figure>
<br/>
```

---
## [References](@id refs)
[^1]:  [F. W. J. Olver and National Institute of Standards and Technology (U.S.), eds., "*NIST Handbook of Mathematical Functions*", Cambridge ; New York: Cambridge University Press : NIST, 2010](https://dlmf.nist.gov/).
[^2]:  J. E. Hansen, ed., *Spherical Near-Field Antenna Measurements*, The Institution of Engineering and Technology, Michael Faraday House, Six Hills Way, Stevenage SG1 2AY, UK: IET, 1988.
[^3]:  J. A. Stratton, "*Electromagnetic Theory*", McGraw-Hill, 1st ed. International series in pure and applied physics, OCLC: 536704, New York, 1941.
[^4]: [T. Limpanuparb, J. Milthorpe: "Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications", arXiv 2014](https://arxiv.org/abs/1410.1748).
