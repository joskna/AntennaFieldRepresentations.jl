# Representation of Radiated Electromangentic Fields via Electromagnetic Surface Currents
The electromagnetic surface equivalence principle, sometimes called Huygens principle, states that the fields ``\bm{E}(\bm{r}) ``, ``\bm{H}(\bm{r}) `` radiated from sources in a volume ``V`` can be alternatively represented by equivalent electric and magnetic surface currents ``\bm{J}_{\mathrm{eq}}(\bm{r}')`` and ``\bm{M}_{\mathrm{eq}}(\bm{r}')`` on a closed surface ``\mathcal{S}`` (with outward pointing normal vector ``\bm{n} ``) enclosing ``V``. 

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
using the dyadic Green's functions ``\bm{\mathcal{G}}_{JE}(\bm{r},\bm{r}')``, ``\bm{\mathcal{G}}_{ME}(\bm{r},\bm{r}')``, ``\bm{\mathcal{G}}_{JH}(\bm{r},\bm{r}')``, and ``\bm{\mathcal{G}}_{MH}(\bm{r},\bm{r}')``. Refer to the [ theory section of dipole representations](@ref "Electromagnetic Fields of Elementary Dipoles") for a definition of the Green's functions.

In a numeric implementation, the surface ``\mathcal{S}`` can be, e.g., approximated by triangles, which serve as a basis for so-called Rao-Wilton-Glisson basis functions.
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