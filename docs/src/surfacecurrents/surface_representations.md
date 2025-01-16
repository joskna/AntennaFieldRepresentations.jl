# Representing Antenna Fields with Equivalent Surface Currents

The [Huygens Principle](@ref Huygens) tells us that the radiated fields of any source distribution confined within a certain volume ``V`` can also be reproduced by equivalent surface currents on the surface ``\partial V`` of the source volume.


In `AntennaFieldRepresentations.jl`, surface current densities are represented by a struct `SurfaceCurrentDensity{P, E, B, C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                                 | Short Description                                                |
| :---------------------------------------- | :--------------------------------------------------------------- |
| `P <: PropagationType`                    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `E <: ElmagType`                          | Can be `Electric` or `Magnetic`                                  | 
| `B <: BEAST.Space`                        | Finite Element space used for expanding the surface currents. Can be any `BEAST.Space{T} where{T <: Real}` from [BEAST.jl](https://github.com/krcools/BEAST.jl/tree/master) .                              |
| `C <: Complex`                            | Element type of the coefficient vector                           |



## Constructors for a `SurfaceCurrentDensity`
To generate a `SurfaceCurrentDensity`, use one of the following constructors:


```julia
SurfaceCurrentDensity{P, E, B, C}(functionspace::B, coefficients::AbstractVector{C}, wavenumber <: Number) where{P <: PropagationType, E <: ElmagType, B <: BEAST.Space{T} where{T <: Real}, C <: Complex}
```
```julia
SurfaceCurrentDensity(::P, ::E, functionspace::B, coefficients::AbstractVector{C}, wavenumber <: Number) where{P <: PropagationType, E <: ElmagType, B <: BEAST.Space{T} where{T <: Real}, C <: Complex}
```
```julia
SurfaceCurrentDensity(::P, ::E, functionspace::B, wavenumber <: Number) where{P <: PropagationType, E <: ElmagType, B <: BEAST.Space{T} where{T <: Real}, C <: Complex}
```
The input arguments for the costructors are

- `P <: PropagationType` : Can be `Radiated`, `Absorbed`, or `Incident`
- `E <: Elmagtype` : Can be `Electric` or `magnetic`
- `coefficients`: Coefficient vector. Can be an `AbstractVector{C}`. If no `coefficients` vector is given, a zero-filled `Vector{C}` of appropriate size will be created for initialization.
- `B <: BEAST.Space` : An object which defines the (triangulated) geometry of the Huygens surface and the space of basis funtions for expanding the surface currents. Refer to the [The `BEAST.Space` Type](@ref beastspace)-section for more detailed information
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## [The `BEAST.Space` Type](@id beastspace)

We don't have to reinvent the wheel. `BEAST.jl` has many relevant surface current expansions already implemented. The documentation for `BEAST.jl` can be found [here](https://krcools.github.io/BEAST.jl/stable/) and it is recommended to have a look at the [tutorial](https://krcools.github.io/BEAST.jl/stable/tutorial/) and the [examples](https://github.com/krcools/BEAST.jl/tree/master/examples) to get used with its general usage.

!!! tip
    It is not strictly necessary to load `BEAST.jl` or `CompScienceMeshes.jl` but the process of generating the geometry and functionspace structs is much simpler once they are loaded via 
    
        using BEAST, CompScienceMeshes
---

Creating a `BEAST.Space`-object contains two steps: Creating the geometry (a mesh) and defining a Finite Element space on this geometry.

### Creating the Geometry (the Mesh)
First, we need to define a mesh. 
The package [`CompScienceMeshes`](https://github.com/krcools/CompScienceMeshes.jl) provides data structures and algorithms for working with simplical meshes in computational science.
We can create a mesh from scratch for very simple meshes, e.g., by
```julia 
using CompScienceMeshes

meshsize = 0.1
sidelengthA = 1.0
sidelengthB = 1.0

Γ = CompScienceMeshes.meshrectangle(sidelengthA, sidelengthB, triangle_sidelength)
```
or 
```julia 
using CompScienceMeshes

radius = 1.0
meshsize = 0.1

Γ = CompScienceMeshes.meshsphere(radius, meshsize)
```


More elaborate meshes can be created with a meshing tool like [`gmsh`](https://gmsh.info/) [^1]. The meshes created by `gmsh` (usually stored in a `.msh`-file) can be loaded for our use with
```julia 
using CompScienceMeshes

Γ = CompScienceMeshes.read_gmsh_mesh("example_box.msh")
```

[^1]: The API of `gmsh` can be accessed via the `julia` package [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl?tab=readme-ov-file). 

### Defining the Finite Element Space

Creating the finite element space is easy in `BEAST.jl`. Several, easy to use constructor functions exist to create many of the most widely used Finite Element spaces for surface current expansion.
To create a set of [Rao-Wilton-Glisson basis functions](@ref rwgbasis) on a previously defined triangular mesh `Γ`, use [^2]
```julia 
using BEAST

β = BEAST.raviartthomas(Γ)
```
---
At the moment, Rao-Wilton-Glisson (RWG) basis functions are the only type of Finite Element space which is tested to work with `AntennaFieldRepresentations.jl`. However, many more Finite Element spaces are defined by `BEAST.jl` which might just "work out of the box" (go ahead and try if you dare and let us know how it worked out for you):

- For rotated RWG (more technically: ``\bm{n} \times ``RWG) functions, use (the `×` symbol can be typed in the REPL by typing `\times` followed by `[TAB]`)
```julia 
using BEAST

β = BEAST.n × BEAST.raviartthomas(Γ)
```
---
- For [Buffa-Christiansen basis functions](https://comptes-rendus.academie-sciences.fr/mathematique/articles/10.1016/j.crma.2004.12.022/) use
 ```julia 
using BEAST

β = BEAST.buffachristiansen(Γ)
```
---
- More Finite Element spaces can be found by consulting the [source code of BEAST.jl](https://github.com/krcools/BEAST.jl/tree/master/src/bases).

[^2]: Rao-Wilton-Glisson basis functions are sometimes called Raviat-Thomas basis functions [(https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions)](https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions).

!!! todo
    Document complete API ...
---