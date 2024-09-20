# [Plane-Wave Representations of Electromagnetic Fields](@id planewavetheory)
This chapter describes how the electromagnetic fields of an antenna can be represented by collections of plane-wave components.

## [Plane-Wave Spectrum Representation of Incident Fields](@id planewave_incident)
In a source-free region of space around the origin incident electormagnetic fields can be represented by an integral of the so-called plane-wave spectrum ``\bm{P}(\vartheta, \varphi)`` over all propagation directions
as
```math
\bm{E}(\bm{r})
=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\bm{P}(\vartheta, \varphi)\,
\mathrm{e}^{-\mathrm{j} \bm{k}(\vartheta, \varphi) \cdot \bm{r}}
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi
```
and
```math
\bm{H}(\bm{r})
=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\dfrac{\bm{k}(\vartheta, \varphi)}{Z_{\mathrm{F}} \lVert\bm{k}(\vartheta, \varphi)\rVert} \times
\bm{P}(\vartheta, \varphi)\,
\mathrm{e}^{-\mathrm{j} \bm{k}(\vartheta, \varphi) \cdot \bm{r}}
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\, ,
```
respectively, where the so-called plane-wave propagation vector direction denotes the direction of propagation
```math
\bm{k}(\vartheta, \varphi)
=
k_0\left(
\sin \vartheta\, \cos \varphi \,\bm{e}_x 
+ 
\sin \vartheta\, \sin \varphi \,\bm{e}_y
+
\cos \vartheta\, \bm{e}_z 
\right)
```
simultaneously encodes the direction of propagation (i.e., the direction of ``\bm{k}``) and the wavelength ``\lambda`` (from its magnitude ``\lVert \bm{k} \rVert =k_0`` via ``\lambda= 2\pi / k_0``) of any plane-wave component.
The plane-wave spectrum ``\bm{P}(\vartheta, \varphi)`` is always perpendicular to the propagation direction ``\bm{k}(\vartheta, \varphi)``.

!!! note
    The above representation of the electric and magnetic field is purely solenoidal, i.e., divergence-free. As such, the representation is only good for representing source-free field solutions of Maxwell's equations. There exists a close relationship between the [spherical vector-wave expansion of incident fields](@ref spherical_incident_expansion) which uses only spherical vector-wave modes ``\bm{F}_{s \ell m}^{(1)}(r, \vartheta, \varphi)`` of incident type and the above representation of the fields as an integral over propagating plane waves.
---

!!! note
    Since the propagation vector ``\bm{k}`` and the corresponding unit vector into the same direction ``\hat{\bm{k}}`` encode the propagation direction of the plane wave, we may use the shorthand notation ``\bm{P}(\hat{\bm{k}})`` to represent the slightly longer expression ``\bm{P}(\vartheta, \varphi)`` whenever convenient (sometimes we may even mix ``\vartheta, \varphi`` with ``\hat{\bm{k}}`` in the same expression). 
    This should not lead to any ambiguities because the relation between ``\hat{\bm{k}}`` and the tuple ``\vartheta, \varphi`` is one-to one, as each encodings uniquely define the same point on the unit sphere.
---

## Far-Field Pattern Representation of Radiated Fields
The radiated field of any current distribution is completely characterized by the radiated far-field pattern 
```math
\bm{F}(\vartheta, \varphi) = \lim \limits_{r\rightarrow \infty} \dfrac{r}{\mathrm{e}^{-\mathrm{j}k_0 r}}\, \bm{E}(r, \vartheta, \varphi)\,,
```
because it is, in prinicple, possible to reconstruct all radiated spherical vector-wave modes from a [spherical expansion of the far fields](@ref spherical_farfield).
Thus, we may represent any radiated field by its far-field pattern.

!!! warning
    Do not confuse the spherical vector-wave modes ``\bm{F}_{s \ell m}^{(1)}(r, \vartheta, \varphi)`` with the far-field pattern ``\bm{F}(\vartheta, \varphi)``!
---

Note, that evaluating the near-fields of a radiated field given in terms of its far-field pattern ``\bm{F}(\vartheta, \varphi)`` is not straightforward. Strictly speaking, it is necessary to reconstruct the spherical vector-wave expansion of the radiated field from its far-field pattern before one can evaluate the spherical expansion at the desired location ``\bm{r}``.


## [Translation of a Far-Field Pattern into a Plane Wave Spectrum in Different Coordinate System](@id planewave_translation)
Fortunately, there is a convenient way to convert the far-field pattern ``\bm{F}(\vartheta, \varphi)`` of a radiated field into an incident plane-wave spectrum ``\bm{P}(\vartheta', \varphi')`` in an observation volume around a new coordinate origin which is translated by the vector ``\bm{R}`` against the original coordinate system.  We use primed coordinates ``r', \vartheta', \varphi' `` and primed vectors ``\bm{k}'`` to distinguish the translated coordinate system from the original coordinate system. The geometrical situation is depicted in the figure below.

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/SphericalTranslation_dark.svg" width="500">
  <source media="(prefers-color-scheme: light)" srcset="../assets/SphericalTranslation.svg" width="500" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic of the geometry for the translation of a far-field representation of radiated fields into a plane-wave spectrum representation of incident fields in a new coordinate system.
  </figcaption>
</figure>
<br/>
```

We have
```math
\bm{P}(\vartheta', \varphi') = T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}}' \cdot \hat{\bm{R}})\, \bm{F}(\vartheta, \varphi)\bigg\rvert_{\vartheta = \vartheta', \varphi = \varphi`}\, ,
```
where 
```math
T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}}' \cdot \hat{\bm{R}})
=
-\sum \limits_{\ell' =0}^{L}\left(-\mathrm{j}\right)^{\ell'} \,
\left(2\ell'+1\right) \, 
\mathrm{h}_{\ell'}^{(2)}(k_0 \lvert \bm{R} \rvert)
\mathrm{P}_\ell'(\hat{\bm{k}}' \cdot \hat{\bm{R}})
```
with ``\hat{\bm{R}}`` the unit vector in the direction of the vector from the original coordinate origin to the new coordinate origin.
Ideally, the number of sumation terms ``L`` should go to ``\infty``, but in practice, we can truncate the summation at a carefully chosen finite ``L``.

The calculated plane-wave spectrum ``\bm{P}(\vartheta', \varphi') `` can be used in the [integral expressions above](@ref planewave_incident) to find the electric and magnetic fields in the observation region.

## [Plane-Wave Representation of the Received Signal of an Antenna](@id planewave_received)
Let ``\hat{\bm{F}}(\vartheta, \varphi)`` be the normalized far field radiated by an antenna such that ``\bm{F}(\vartheta, \varphi)=  a\,\hat{\bm{F}}(\vartheta, \varphi)`` is the actual radiated far-field when its transmit port is excited by a signal with wave amplitude ``a \in \mathbb{C} \, \sqrt{\mathrm{W}}``.
To determine the received signal ``b \in \mathbb{C} \, \sqrt{\mathrm{W}}`` of this antenna in receive mode under a certain illumination of an incident field, it is convenient to represent the incident field via its plane-wave spectrum ``\bm{P}(\vartheta, \varphi)``. In this case, the received signal is given by (remember that ``\hat{\bm{k}}`` depends on ``\vartheta, \varphi``)
```math
b=
\int \limits_{0}^{2\pi}
\int \limits_{0}^{\pi}
\bm{P}(\hat{\bm{k}})\cdot
\hat{\bm{F}}\left(- \hat{\bm{k}} \right)
\, \sin \vartheta\,
\mathrm{d}\vartheta\, \mathrm{d}\varphi\, .
```

Combining this result with the findings from the [section above](@ref planewave_translation), we can express the ``S_{21}`` parameter measured for the transmission between two antennas as
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
where ``\hat{\bm{F}}_1\left(\hat{\bm{k}} \right)`` is the normalized far-field pattern of antenna 1 (the transmit antenna), ``\hat{\bm{F}}_2\left( - \hat{\bm{k}} \right)`` is the inverted far-field pattern of antenna 2 (the receive antenna) and ``T_L(k_0 \lvert \bm{R} \rvert, \hat{\bm{k}} \cdot \hat{\bm{R}})`` is the translation operator known from [previous sections](@ref planewave_translation) with ``\bm{R}`` denoting the vector separating the centers of the two antennas.
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

