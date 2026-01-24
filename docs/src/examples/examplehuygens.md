
## Equivalent Surface Currents on a Huygens Surface Around a Huygens Radiator

In this example, we will define equivalent currents on a surface (the so-called Huygens surface) enclosing a distribution of sources such that the equivalent currents will radiate the same field as the original sources. 
```@contents
Pages = ["examplehuygens.md"]
Depth = 3
```

### Setup of the Original Sources
For the original sources, we take a suitable combination of a `z`-oriented electric dipole and a `y`-oriented magnetic dipole to generate a source distribution which predominantly radiates into the positive `x`-direction: 

```jldoctest huygensexample ; output=false
using AntennaFieldRepresentations
using LinearAlgebra

λ = 2.0 # set an arbitrary wavelength
k0 = λ / (2π) # determine the corresponding wavenumber

# setup a "collection" of electrical dipoles (= Hertzian dipoles) consisting of only one element
# oriented in z-direction at R= [0.5, 0.5, 0.5]] with unit excitation
eldips = HertzArray([[0.5, 0.5, 0.5]], [complex.([0.0, 0.0, 1.0])], [complex(1.0)], k0)

# setup a "collection" of corresponding magnetic dipoles (= Fitzgerald dipoles) consisting of only one element
# oriented in negative y-direction at R= [0.5, 0.5, 0.5]] with excitation scaled by the free space impedance Z₀ 
magdips = FitzgeraldArray(
    [[0.5, 0.5, 0.5]],
    [complex.([0.0, -1.0, 0.0])],
    AntennaFieldRepresentations.Z₀ * [complex(1.0)],
    k0,
)

length(magdips)
# output

1

```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/HuygensDipoles.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/HuygensDipoles.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Original field representation: A combination of a Hertzian (green) and a Fitzgerald dipole (purple). 
  </figcaption>
</figure>
<br/>
```

Such a combination of an electric and a magnetic dipole is sometimes called a Hygens radiator.
We can confirm that the Huygens radiator mainly radiates in `x`-direction by converting the antenna field representation into a far-field representation:
```jldoctest huygensexample ; output=false
# convert each DipoleArray into a plane wave representation
pwe_el = changerepresentation(PlaneWaveExpansion, eldips)
pwe_mag = changerepresentation(PlaneWaveExpansion, magdips)

# copy the structure of the plane wave representation
pwe_dipoles = copy(pwe_el)

# the plane-wave coefficients of the total field are the sum of the plane wave coefficients of the separate dipole fields.
pwe_dipoles .= pwe_el + pwe_mag
 #evaluates to 106.4275744120224
round(norm(pwe_dipoles), digits =2)

# output

106.43

```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/HuygensFarfield.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/Huygensfarfield.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The far field of the Huygens radiator has its main beam in x-direction. 
  </figcaption>
</figure>
<br/>
```

### Setup of the Equivalent Surface Current Density Model

Now, we define a triangulated Huygens surface with a suitable set of basis functions for the equivalent surface currents.
The creation of the mesh is simple: Define an array of 3D arrays to define the vertices of the mesh and define a second array of integer triples which define the three vertices which generate a triangular face of the mesh like so:
```jldoctest huygensexample ; output=false
vertices = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 1.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [0.5, 0.5, 1.0],
    [0.5, 1.0, 0.5],
    [1.0, 0.5, 0.5],
]

faces = [
    # front
    [1, 5, 10],
    [5, 6, 10],
    [6, 2, 10],
    [2, 1, 10],

    # bottom
    [1, 9, 5],
    [5, 9, 7],
    [7, 9, 3],
    [3, 9, 1],

    # left
    [1, 2, 11],
    [2, 4, 11],
    [4, 3, 11],
    [3, 1, 11],

    # right
    [5, 7, 14],
    [7, 8, 14],
    [8, 6, 14],
    [6, 5, 14],

    # back 
    [7, 13, 8],
    [8, 13, 4],
    [4, 13, 3],
    [3, 13, 7],

    # top
    [2, 6, 12],
    [6, 8, 12],
    [8, 4, 12],
    [4, 2, 12],
]

# output

24-element Vector{Vector{Int64}}:
 [1, 5, 10]
 [5, 6, 10]
 [6, 2, 10]
 [2, 1, 10]
 [1, 9, 5]
 [5, 9, 7]
 [7, 9, 3]
 [3, 9, 1]
 [1, 2, 11]
 [2, 4, 11]
 ⋮
 [6, 5, 14]
 [7, 13, 8]
 [8, 13, 4]
 [4, 13, 3]
 [3, 13, 7]
 [2, 6, 12]
 [6, 8, 12]
 [8, 4, 12]
 [4, 2, 12]
```
The rest of the heavy lifting for the generation of basis functions is offloaded to `CompScienceMeshes.jl` and `BEAST.jl`:
```jldoctest huygensexample ; output=false
using CompScienceMeshes: Mesh
using BEAST: raviartthomas
using StaticArrays: SVector


mesh = Mesh(SVector{3}.(vertices), SVector{3}.(faces))
Γ = raviartthomas(mesh)


excitations = zeros(ComplexF64, length(Γ))

# output

36-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
```
Finally, we transform the set of RWG basis functions (aka Raviart-Thomas basis functions) into a `SurfaceCurrentDensity`:
```jldoctest huygensexample ; output=false
currents = SurfaceCurrentDensity{Radiated,Electric,typeof(Γ),ComplexF64}(Γ, excitations, k0)

# output

36-element SurfaceCurrentDensity{Radiated, Electric, BEAST.RTBasis{Float64, Mesh{3, 3, Float64}, SVector{3, Float64}}, ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
```

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/HuygensCurrentsEmpty.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/HuygensCurrentsEmpty.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The basis functions for the equivalent surface currents are defined as RWG functions on a triangulated, box shaped surface around the original soures.
  </figcaption>
</figure>
<br/>
```

### Finding the Correct Excitation Coefficients for the Equivalent Surface Current Density

We are half way done with replacing the original sources by an equivalent electrical surface current density. 
We are missing the correct excitation coefficients for the RWG-basis functions on the surface (at the moment, all excitation coefficients are zero).
To find the correct coefficients, we have to ensure that the radiated fields of the `SurfaceCurrentDensity` are (almost) identical to the radiated fields of the original sources.
Fortunately, we only have to enforce the radiated far-fields to be identical. The radiated near-fields then must also coincide[^1].

[^1]: In principle, evanescent fields can lead to deviations between the very-near fields of the equivalent `SurfaceCurrentDensity` and the original sources. For the sake of simplicity we neglect all technical issues with evanescent fields.

To this end, we store the radiated far fields of each individual RWG-basis function in a matrix `A_farfield` (i.e., every coloumn of the matrix `A_farfield` corresponds to the radiated far-field of a single RWG-basis function)[^2]:
```jldoctest huygensexample ; output=false , filter = r"(.*\R*.*)*" => s""
# initialize zero-filled matrix
A_farfield = zeros(ComplexF64, length(pwe_dipoles), length(currents))

# fill the matrix
for k in 1:length(currents)
    currents .= 0
    currents[k] = 1

    A_farfield[:, k] .= changerepresentation(PlaneWaveExpansion, currents; samplingstrategy=pwe_dipoles.samplingstrategy)
end


# output

```
[^2]: This is an illustrative example. For realistic cases, the matrix size usually prohibits to calculate all matrix elements explicitly. In such scenarios, we can use `AntennaFieldRepresentations.jl`'s  `TransmitMap`s which neatly integrate with iterative solvers, e.g., from [`IterativeSolvers.jl`](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl).


Finally, we solve for the desired excitation coefficients by solving a linear system of equations. The desired right-hand side is the vector of far-field coefficients stored in the variable `pwe_dipoles` (remember that all `AntennaFieldRepresentations` behave like an `AbstractVector`).

```jldoctest huygensexample ; output=false 
currents .= A_farfield \ pwe_dipoles
length(currents)

# output

36
```


```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/HuygensCurrents.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/HuygensCurrents.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Visualization of the equivalent electric surface currents on the Huygens box which radiates the same field as the original sources.
  </figcaption>
</figure>
<br/>
```

The newly found equivalent surface currents provide a complete representation of the radiated fields and can be utilized for various post-processing applications. As a final step for this example, we want to verify that the corresponding far fields actually match the far fields of the original sources:

```jldoctest huygensexample ; output=false, filter = r"(\p{L}*\R*.)*" => s""
using LinearAlgebra
pwe_currents = changerepresentation(PlaneWaveExpansion, currents; samplingstrategy=pwe_dipoles.samplingstrategy)

norm(pwe_dipoles - pwe_currents)/norm(pwe_dipoles) < 2e-6 # true

# output

```