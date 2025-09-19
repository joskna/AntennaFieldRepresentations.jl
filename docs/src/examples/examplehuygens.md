
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

# output

1-element FitzgeraldArray{Float64, ComplexF64}:
 376.73031366686996 + 0.0im

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

# output

196-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
      -2.507461811441103 + 12.29508945871003im
      -3.554309025967161 + 15.543864561707764im
     -3.8082481186806953 + 17.86271487494979im
     -3.0247253020703613 + 18.844171281879387im
     -1.4753417863809353 + 18.20446931526832im
     0.17919273579267136 + 15.944049928679517im
      1.2643276978190043 + 12.484312756657681im
     -2.5091139861471126 + 11.328608703711998im
      -3.864245379246303 + 14.493743842503624im
      -4.431616542349483 + 16.74255276176831im
                         ⋮
     -0.7442018868542088 - 5.482115441164023im
     -1.1210998092957607 - 6.991769960174524im
     -0.6821359032386382 + 3.870046675223574im
     -0.5130571603587533 + 3.027078197822951im
    -0.22205406878586972 + 1.6656322337757892im
 -1.8829602583702385e-17 + 2.528272837396797e-16im
    0.005627840829237915 - 1.68035921012886im
      -0.209054152967389 - 3.063123638813411im
       -0.50023468457528 - 3.8977347163761804im

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

# output

36-element SurfaceCurrentDensity{Radiated, Electric, BEAST.RTBasis{Float64, Mesh{3, 3, Float64}, SVector{3, Float64}}, ComplexF64}:
 -0.12599404604722858 - 0.0014105146213536307im
 7.492148395128541e-9 - 1.5757718018659355im
   0.1259940310048704 + 1.5757717167798422im
  -0.2523956434097737 - 3.158767714407309im
 9.056014724688309e-9 - 4.990083599449007e-9im
   0.2523956254112038 + 4.943777054914383e-9im
 7.507966032954928e-9 + 1.5757717936312878im
  -0.1259940460643345 + 0.0014105229089590705im
  9.08089751949105e-9 - 4.997434211502562e-9im
  -0.1259940309773886 + 1.5757717084946292im
                      ⋮
  -0.2523956252740142 + 3.158767714302966im
  0.12599404610205847 - 1.575771716818354im
  -0.1259940309474584 + 0.0014105145865028193im
 7.606307412308501e-9 - 1.5757718019492395im
  -0.2523956435497995 + 5.040879940648099e-9im
  0.12599403094719086 + 0.0014105229383432995im
 -0.12599404608976275 - 1.575771708469698im
  -0.2523956253018469 - 3.158767724377384im
 7.605682896331188e-9 + 1.5757717935542572im

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