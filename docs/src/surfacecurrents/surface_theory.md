# [Huygens Principle and Equivalent Electromagnetic Surface Currents](@id Huygens)
The [electromagnetic surface equivalence principle](https://en.wikipedia.org/wiki/Surface_equivalence_principle), sometimes called Huygens principle, states that the fields ``\bm{E}(\bm{r}) ``, ``\bm{H}(\bm{r}) `` radiated from sources in a volume ``V`` can be alternatively represented by equivalent electric and magnetic surface currents ``\bm{J}_{\mathrm{eq}}(\bm{r}')`` and ``\bm{M}_{\mathrm{eq}}(\bm{r}')`` on a closed surface ``\mathcal{S}`` (with outward pointing normal vector ``\bm{n} ``) enclosing ``V\,.`` 

The corresponding fields fields can be calculated at any point ``\bm{r} `` outside of ``\mathcal{S} `` from the equivalent surface currents via 
```math
\bm{E}(\bm{r})=
\iint \limits_{\mathcal{S}}
\left[
\bm{\mathcal{G}}_{JE}
(\bm{r},\bm{r}')
\cdot
\bm{J}_{\mathrm{eq}}(\bm{r}')
+
\bm{\mathcal{G}}_{ME}
(\bm{r},\bm{r}')
\cdot
\bm{M}_{\mathrm{eq}}(\bm{r}')
\right]\, \mathrm{d}a' 
```
and
```math
\bm{H}(\bm{r})=
\iint \limits_{\mathcal{S}}
\left[
\bm{\mathcal{G}}_{JH}
(\bm{r},\bm{r}')
\cdot
\bm{J}_{\mathrm{eq}}(\bm{r}')
+
\bm{\mathcal{G}}_{MH}
(\bm{r},\bm{r}')
\cdot
\bm{M}_{\mathrm{eq}}(\bm{r}')
\right]\, \mathrm{d}a' 
```
using the dyadic Green's functions 
```math
\bm{\mathcal{G}}_{JE}=
-\mathrm{j} k_0\, Z_\mathrm{F}
\left[
\left(\mathbf{I} +\dfrac{1}{k^2} \nabla \nabla\right) g_0(\bm{r},\bm{r}')
\right]\,,
```
```math
\bm{\mathcal{G}}_{ME}=
-\nabla g_0(\bm{r},\bm{r}') \times \mathbf{I}\,,
```
```math
\bm{\mathcal{G}}_{JH}=
\nabla g_0(\bm{r},\bm{r}') \times \mathbf{I}\,,
```
and
```math
\bm{\mathcal{G}}_{MH}=
-\mathrm{j} \dfrac{k_0}{Z_\mathrm{F}}
\left[
\left(\mathbf{I} +\dfrac{1}{k^2} \nabla \nabla\right) g_0(\bm{r},\bm{r}')
\right]
```
which are derived from the scalar greens function
```math
 g_0(\bm{r},\bm{r}')
=
\dfrac{\mathrm{e}^{-\mathrm{j}k_0\, ||\bm{r}-\bm{r}'||}}{4\pi\, ||\bm{r}-\bm{r}'||}\, .
```

## [Rao-Wilton-Glisson Basis Functions](@id rwgbasis)
In a numeric implementation, the surface ``\mathcal{S}`` can be, e.g., approximated by triangles.
On these triangles, we can define a set of basis functions ``\bm{\beta}_n(\bm{r})`` such that the equivalent electric surface current density [^1] ``\bm{J}_{\mathrm{eq}}(\bm{r})`` on the triangulated surface can be expressed as a linear superposition
```math
\bm{J}_{\mathrm{eq}}(\bm{r}) = \sum \limits_{n=1}^N
j_n\, \bm{\beta}_n(\bm{r})
```
with some coefficients ``j_n \in \mathbb{C}\,  \mathrm{A}/\mathrm{m}``.

A very common choice of basis functions on a triangulated surface are the so-called Rao-Wilton-Glisson (RWG) basis functions [^2].
RWG basis functions are sometimes called [Raviart-Thomas basis functions](https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions).
They are defined on each pair of adjacient triangles which share a common edge. The current density for an RWG-basis function crosses the common edge and at each point ``\bm{r}`` inside one of the two triangles, the current density points away from or towards the vertex of the triangle which is not part of the common edge (the vertex is denoted by ``\bm{r}_{\mathrm{op}-}`` or ``\bm{r}_{\mathrm{op}+}`` depending of whether the current density points away from or towards this vertex).

Formally, the surface current density of the ``n``th RWG-basis function is given by
```math
\bm{\beta}_n(\bm{r})
=\pm
\dfrac{\bm{r}_{\mathrm{op}\pm} -\bm{r}}{2A_\Delta}\, ,
```
where ``A_\Delta`` denotes the surface area of the triangle in which ``\bm{r}`` is located. A schematic representation of an RWG function is depicted below.

[^1]: The same procedure also works for the magnetic surface current denisity.
[^2]: Rao, Sadasiva, Donald Wilton, and Allen Glisson. "Electromagnetic scattering by surfaces of arbitrary shape." IEEE Transactions on antennas and propagation 30.3 (1982): 409-418. [doi: 10.1109/TAP.1982.1142818](https://ieeexplore.ieee.org/document/1142818)

```@raw html
<figure>
<picture>
  <source  srcset="../assets/RWG.png" width="200">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Schematic representation of Rao-Wilton Glisson function.
  </figcaption>
</figure>
<br/>
```

The total surface current density is created by superimposing weighted RWG-basis functions for each triangle pair, such that a current density is defined at each point on the triangulated surface as in the figure below.

```@raw html
<figure>
<picture>
  <source  srcset="../assets/SurfacePlot.png" width="600">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Electric surface currents represented by Rao-Wilton-Glisson functions on spherical surface.
  </figcaption>
</figure>
<br/>
```

With the known surface current density on the complete triangulated surface, we can use numeric integration techniques to (approximately) evaluate the integrals for calculating the electic and magnetic fields.


## Non-Uniqueness of Equivalent Surface Current Representations

The representation of antenna fields by equivalent surface current densities is ambiguous in the sense that various different surface currents can generate exactly the same or approximately the same fields (in the volum outside of the original source volume). 

The linearity of Maxwell's equations implies the existence of non-radiating currents or weakly radiating surface currents [^3] in this case. If two different surface current distributions generate the same field, their difference (i.e., a certain other current distribution) must result in a zero exterior field.

!!! note
    Non-radiating surface currents can be added to any surface current density without affecting the exterior fields.
---

[^3]: By weakly radiating currents, we understand a surface current density which does not contribute significantly to the radiated far-fields despite the magnitude of the currents being in the same order of magnitude to other currents which do contribute to the far field pattern.

### Non-Radiating Currents
Non-radiating surface currents have zero exterior fields but non-zero interior fields. 
More technically, let ``V`` be the interior (original source-) volume and ``\partial V`` its boundary with an outward pointing unit normal vector ``\bm{n}``.

!!! todo
    Complete section about surface equivalence principle and non-radiating, weakly radiating and Love currents
### Weakly Radiating Currents

### Love Currents