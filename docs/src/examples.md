# Examples

!!! warning
    This chapter refers to an old state of the package and needs rewriting!! 
---

## Convert Radiated Fields of Fitzgerald Dipoles Into Other Representations
First, let us investigate the spherical mode expansion of fields which are radiated by a collection of Fitzgerald dipoles. 

`AntennaFieldRepresentations` is used to convert the dipole representation into a radiating spherical vector wave expansion. 
TODO: add example

The order of the mode expansion (i.e., th number of considered spherical mode coefficients to represent the radiated fields) is based on the largest separation of any dipole to the coordinate origin. The default value is chosen very conservatively, leading to many close-to-zero coefficients as can be seen by visualizing the corresponding spherical mode expansion by
TODO: add example


```@raw html
<figure>
<picture>
  <source  srcset="../assets/sphModes.svg" width="800">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Mode spectrum for the radiated field of the collection of Fitzgerald dipoles.
  </figcaption>
</figure>
<br/>
```
A larger expansion order (i.e., the maximum occuring mode order ℓ) usually leads to more accurate representations at the cost of additional storage and computation requirements. The expansion order can be passed as an additional argument to the `convertrepresentation` function. In the following example, the expansion order is limited to an expansion order of ℓ≤24.
TODO: add example

```@raw html
<figure>
<picture>
  <source  srcset="../assets/sphModes_reduced.svg" width="800">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Reduced mode spectrum for the radiated field of the collection of Fitzgerald dipoles.
  </figcaption>
</figure>
<br/>
```


Next, we can calculate the radiated far fields and store the far-field pattern in a `FarFieldPattern.`
TODO: add example


For a higher resolution, the mode order for the far-field representation can be increased.
TODO: add example

```@raw html
<figure>
<picture>
  <source  srcset="../assets/ff_dipoles_highres.png" width="800">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Far-field pattern of the collection of Fitzgerald dipoles.
  </figcaption>
</figure>
<br/>
```

When the far-field pattern is computed from a spherical wave expansion, the mode order for the far-field representation is chosen identical to the mode order of the spherical wave expansion.
TODO: add example

```@raw html
<figure>
<picture>
  <source  srcset="../assets/ff_sph.png" width="800">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Far-field pattern of the collection of Fitzgerald dipoles.
  </figcaption>
</figure>
<br/>
```

