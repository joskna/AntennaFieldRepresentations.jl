# Conversions Between Field Representations
Field representations can be converted into one another using the `convertrepresentation` function. 
The fields which are evaluated at any given location from the field representation before and after the conversion should be approximately the same if the location is in a region which is valid for both representations.
This section is intended to clarify which oenversions are possible and to identify the common regions of validity for the converted field representations. 
 

!!! note
    In general, a representation of an _incident_ field can only be converted into another _incident_ field representation and a _radiating_ representation can only be converted into another _radiating_ representation.

    Use the `translate`- methods to convert _radiating_ representations into their corresponding _incident_ representation in a translated coordinate system.
---

!!! note
    The transmission, e.g. between two antennas, should be exactly reciprocal. Therefore, in theory, the order of arguments in a `transmission`-command should be exchangable without affecting the result. The calculation methods differ, however. Internally, the fields of the source representation (first argument) are always converted into a representation suitable to be received by the receiver representation (second argument), possibly breaking the exact reciprocity.
---
!!! warning
    Different field representations may have different regions of validity. The result of `convertrepresentation` may only approximately evaluate to the same electric or magnetic fields as the original representation in a region of space which is valid for the old *and* the new representation.
---


## Collection of Dipoles ``\rightarrow`` Radiating Spherical Wave Expansion
Radiated fields can be represented by a vector of dipoles `dipoles::Array{<:AbstractDipole,1}` (the vector may contain any number of Hertz dipoles or Fitzgerald dipoles). 
In theory, the fields of the dipoles can be evaluated at any location which does not exactly coincide with one of the dipole locations. Exactly at the dipole locations, the fields have a singularity. In practice, if the fields shall be evaluated too close to a dipole location (much less than ``\lambda/1000`` distance) overflows may happen.

The field representation can be converted into a radiating spherical wave expansion by
TODO: add example

In the above example, the coefficients of the new spherical wave expansion are represented by complex numbers of type `ComplexF64`. If single precision is sufficient, one may use
TODO: add example


A more effective way to control the accuracy of the spherical expansion is to control the expansion order, i.e., the number of mode coefficients which is used for the field representation. The expansion order of the resulting mode expansion is controlled by the Lmax argument
TODO: add example

Also have a look at [Convert Radiated Fields of Fitzgerald Dipoles Into Other Representations](@ref) in the [Examples](@ref) section of this documentation.

In theory, the radiating spherical wave expansion is valid outside a sphere centered at the origin with a radius large enough to enclose all source dipoles. In practice, numerical instabilities prevent a reliable evaluation of the fields at locations ``\bm{r}`` with ``||k_0 \bm{r}|| \ll \ell_{\mathrm{max}}``, where ``\ell_{\mathrm{max}}`` denotes the largest mode order ``\ell`` contained in the radiating spherical mode expansion.

!!! note
    We can directly convert dipole based representations into other field representations, but not the other way around.

    (Equivalent) current based field representations, such as a collection of Hertz or Fitzgerald dipoles, are not unique: Different current distributions can lead to the exact same radiated fields. Since there is no unique current distribution for a given field,there is no simple conversion from a field representation to a dipole based representation.
---
## Collection of Dipoles ``\rightarrow`` Far-Field Pattern

## Collection of Dipoles ``\rightarrow`` Incident Spherical Wave Expansion 

## Radiating Spherical Wave Expansion ``\rightarrow`` Far-Field Pattern
## Far-Field Pattern ``\rightarrow`` Radiating Spherical Wave Expansion

## Incident Spherical Wave Expansion ``\rightarrow`` Plane Wave Spectrum
## Plane Wave Spectrum ``\rightarrow`` Incident Spherical Wave Expansion

