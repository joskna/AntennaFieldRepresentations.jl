# Field Representations in Transformed Coordinate Systems

In many applications, it is useful to rotate or translate an `AntennaFieldRepresentation` into a new position. Of course, this is equivalent to a `AntennaFieldRepresentation` represented in a inversely rotated or translated coordinate system.


## Rotated Coordinate System

One can use the [`rotate`](@ref) function or its in-place [`rotate!`](@ref) version to rotate an `AntennaFieldRepresentation` into its new position.
The `rotate(aut_field :: AntennaFieldRepresentation, χ::Number, θ::Number, ϕ::Number)` function rotates the field representation where the rotation is defined by the Euler angles `χ`, `θ`, `ϕ`.
This is equivalent to the original representation being represented in a rotated coordinate system, where the coordinate axes of the original coordinate system must be rotated around the Euler angles -`χ`, -`θ`, -`ϕ` to get the rotated coordinate frame.

The rotations are performed according to the z,y,z- sequence defined by extrinsic rotations. Extrinsic rotations are elemental rotations that occur about the axes of the fixed coordinate, i.e., not the intrinsic coordinates of the rotated object. First, the object is rotated about the global     
z-axis by `χ`. Then, the object is rotated about the global y-axis by `θ`. Finally, the object is rotated about the global z-axis by `ϕ`. 

Let's observe the rotations in action in the following example.
Consider a farfield pattern as shown in the figure below (expand "Setup Code" to see how the input data is generated):

```@raw html
<details closed><summary>Setup Code</summary>
```

```jldoctest rotateexamples ; output=false
#########################################
#                   Setup
#                     |
#                     V
#########################################
using AntennaFieldRepresentations

function generate_AUTdips(xvec::Array{Float64,1}, yvec::Array{Float64,1}, zvec::Array{Float64,1}, k0::Float64)

    nx = length(xvec)
    ny = length(yvec)
    nz = length(zvec)

    ycenter = (maximum(yvec) + minimum(yvec)) / 2
    ysize = maximum([abs(maximum(yvec) - ycenter), abs(minimum(yvec) - ycenter)])

    ndips = nx * ny * nz
    positions = Vector{Vector{Float64}}(undef, ndips)
    magnitudes = Vector{ComplexF64}(undef, ndips)

    # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
                # println(index)
                positions[index] = [xvec[kkk], yvec[kk], zvec[k]]
                magnitudes[index] = mag
                # println(pos)


            end
        end
    end
    return HertzArray(positions, [complex.([0.0, 0.0, 1.0]) for k = 1:length(positions)], magnitudes, k0)

end


using LinearAlgebra


Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = generate_AUTdips(collect(-0.25λ:λ/4:0λ), collect(-0.5λ:λ/4:0.5λ), collect(-1λ:λ/4:1λ), k0)

pwe = rotate(changerepresentation(PlaneWaveExpansion, dipoles), 0.0, pi / 2, 0.0)
round.(pwe, digits=2)
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

4096-element Vector{ComplexF64}:
   133.27 - 0.28im
   565.85 - 6.4im
   901.33 - 24.95im
   741.71 - 37.95im
    11.31 - 0.92im
   -993.1 + 118.04im
  -1812.5 + 295.13im
 -2109.39 + 450.36im
 -1826.43 + 494.96im
 -1154.86 + 387.6im
          ⋮
  -272.69 + 1006.23im
  -281.75 + 1319.64im
  -220.95 + 1356.97im
  -111.43 + 937.48im
     3.78 - 46.33im
    79.22 - 1548.24im
    92.13 - 3328.13im
    56.38 - 4989.23im
     13.1 - 6097.33im

``` 


```@raw html
</details>
```

We start with a farfield pattern `pwe` represented as `PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}` which points into the negative z-direction:
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/pwe.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/pwe.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Original far field pattern pointing into negative z-direction. 
  </figcaption>
</figure>
<br/>
```


Now, let us rotate this pattern by π / 4 along θ:
```jldoctest rotateexamples ; output=false
θ = π / 4

pwe_rotθ = rotate(pwe, 0.0, θ, 0.0)
round.(pwe_rotθ, digits = 2)
# output

4096-element Vector{ComplexF64}:
 -2065.64 + 398.34im
 -1546.83 + 224.55im
  -606.62 + 63.01im
   344.58 - 23.84im
   874.82 - 36.22im
   814.63 - 16.77im
    384.2 - 2.65im
    33.03 + 0.01im
    90.64 - 0.12im
    504.0 - 4.88im
          ⋮
   239.97 - 219.07im
   243.24 - 258.38im
   225.85 - 279.39im
   187.02 - 270.23im
   126.47 - 214.53im
    46.98 - 94.27im
   -42.81 + 102.68im
  -127.21 + 369.65im
  -185.91 + 665.48im
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/pwe_rottheta.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/pwe_rottheta.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Rotated far-field pattern by  θ = π / 4 around the y-axis. 
  </figcaption>
</figure>
<br/>
```


Now, let us rotate this pattern by π / 4 along θ and also by π / 3 around ϕ:
```jldoctest rotateexamples ; output=false
θ = π / 4
ϕ = π / 3

pwe_rotθ = rotate(pwe, 0.0, θ, ϕ)
length(pwe_rotθ)

# output

4096
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/pwe_rotthetaphi.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/pwe_rotthetaphi.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Rotated far-field pattern by  θ = π / 4 around the y-axis and then by ϕ = π / 3 around the z-axis. 
  </figcaption>
</figure>
<br/>
```

Finally, observe what happens when the far-field pattern is initially rotated by χ = π / 2 around the z-axis before performing the above sequence of rotations.
The far-field pattern is now rotated around its main lobe (which was originally pointing into the negative z-direction).
This is very useful for rotating the probe polarization in spherical measurement setups.
```jldoctest rotateexamples ; output=false
χ = π / 2
θ = π / 4
ϕ = π / 3

pwe_rotχθϕ = rotate(pwe, χ, θ, ϕ)
length(pwe_rotχθϕ)

# output

4096
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/pwe_chirotthetaphi.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/pwe_chirotthetaphi.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
     Rotated far-field pattern by χ = π / 2 around the z-axis, then by  θ = π / 4 around the y-axis, and then by ϕ = π / 3 around the z-axis. 
  </figcaption>
</figure>
<br/>
```

Since the rotation of a `PlaneWaveRepresentation` is ultimately based on interpolating the stored samples to new sampling locations after the rotation, there are additional parameters to control the accuracy of the interpolation.
Only for the special case of `aut_field` being a `PlaneWaverepresentation`, the user can specify the iterpolation orders along `θ` and `ϕ` by using the optional keyword arguments `orderθ` and `orderϕ` in the function `rotate(aut_field :: PlaneWaveRepresentation, χ::Number, θ::Number, ϕ::Number; orderθ=12, orderϕ=12)`. See the following example: 
```jldoctest rotateexamples ; output=false
χ = π / 2
θ = π / 4
ϕ = π / 3

# Default is orderθ=12, orderϕ=12
pwe_rotχθϕ_quick = rotate(pwe, χ, θ, ϕ; orderθ=8, orderϕ=8)
length(pwe_rotχθϕ_quick)

# output

4096
```

In the end, the user has to find the tradeoff between a faster or a more accurate calculation
```jldoctest rotateexamples ; output=false
using LinearAlgebra

# Evaluates to 7.537201623735068e-6 :
round(norm(pwe_rotχθϕ_quick - pwe_rotχθϕ) / norm(pwe_rotχθϕ), digits=8)

# output

7.54e-6
```

### `RotateMap`s
If the same rotation will be applied multiple times to the same type of `AntennaFieldRepresentation`, it is beneficial to define a corresponding `RotateMap`. A `RotateMap` represents the linear operator which takes the coefficients of an `AntennaFieldRepresentation` as input and returns the coefficients of the corresponding rotated `AntennaFieldRepresentation`.  

Analoguous to the `rotate` command, a `RotateMap` is constructed via the constructor[^1]

```julia
R = RotateMap(aut_field::AntennaFieldRepresentation, χ::Number, θ::Number, ϕ::Number)
```
[^1]: Only for the special case of `aut_field` being a `PlaneWaverepresentation`, the user can specify the iterpolation orders along `θ` and `ϕ` by using the optional keyword arguments `orderθ` and `orderϕ` in the constructor `RotateMap(aut_field :: PlaneWaveRepresentation, χ::Number, θ::Number, ϕ::Number; orderθ=12, orderϕ=12)`.

All `RotateMap`s are a subtype of the abstract type [`OperationMap`](@ref operationmaps_linmap), i.e., they behave as liniear maps.
```julia
R = RotateMap(pwe, χ, θ, ϕ)

rotated_pwe = R * pwe
```

The transpose, adjoint, and inverse operators are obtained by applying the `transpose`, `adjoint`, or `inverse` command, respectively.
```julia
Rᵀ = transpose(R)

Rᴴ = adjoint(R)

R⁻¹ = inverse(R)
retrieved_pwe = R⁻¹ * rotated_pwe

norm(retrieved_pwe - pwe) / norm(pwe)
```


## Spatially Shifted Coordinate System

One can use the [`spatialshift(aut_field::AntennaFieldRepresentation, R::AbstractVector)`](@ref) function or its in-place [`spatialshift!`](@ref) version to move (or "shift in space") an `AntennaFieldRepresentation` along the vector `R` to its new location.

Notice that the terms "translation" or "translate" has been avoided in `AntennaFieldRepresentations.jl` due to potential confusion with the commonly so-called "translation" operator in the MLFMM. Instead, the term "spatial shift" is used to express the action of moving an `AntennaFieldRepresentation` along a vector `R` (i.e., in particular the `PropagationType` of the `AntennaFieldRepresentation` is not affected by a spatial shift).

### `SpatialShiftMap`s
If the same spatial shift will be applied multiple times to the same type of `AntennaFieldRepresentation`, it is beneficial to define a corresponding `SpatialShiftMap`. A `SpatialShiftMap` represents the linear operator which takes the coefficients of an `AntennaFieldRepresentation` as input and returns the coefficients of the corresponding spatially shifted `AntennaFieldRepresentation`.  

Analoguous to the `spatialshift` command, a `SpatialShiftMap` is constructed via the constructor
```julia
T = SpatialShiftMap(aut_field::AntennaFieldRepresentation, R::AbstractVector)
```

All `SpatialShiftMap`s are a subtype of the abstract type [`OperationMap`](@ref operationmaps_linmap), i.e., they behave as liniear maps.

The transpose, adjoint, and inverse operators are obtained by applying the `transpose`, `adjoint`, or `inverse` command, respectively.

## Transfer of Radiated Representations into Incident Representations
Radiated spherical mode expansions and radiated plane-wave expansions (represented with respect to a certain coordinate origin) can be expressed as incident expansions represented in a different coordinate system which is shifted by a vector `R` from the original coordinate origin. This process is commonly called "translation" but in order to avoid confusion with a "spatial shift", the process is called "transfer" in `AntennaFieldRepresentations`.

 One can use the `transfer(aut_field::AntennaFieldRepresentation{Radiated}, R::Abstractvector)` method or its in-place version `transfer` to perform the transfer of a `Radiated` field expansion into an `Incident` one (expand "Setup Code" to see how the input data is generated):

```@raw html
<details closed><summary>Setup Code</summary>
```

```jldoctest transferexamples ; output=false
#########################################
#                   Setup
#                     |
#                     V
#########################################
using AntennaFieldRepresentations

function generate_AUTdips(xvec::Array{Float64,1}, yvec::Array{Float64,1}, zvec::Array{Float64,1}, k0::Float64)

    nx = length(xvec)
    ny = length(yvec)
    nz = length(zvec)

    ycenter = (maximum(yvec) + minimum(yvec)) / 2
    ysize = maximum([abs(maximum(yvec) - ycenter), abs(minimum(yvec) - ycenter)])

    ndips = nx * ny * nz
    positions = Vector{Vector{Float64}}(undef, ndips)
    magnitudes = Vector{ComplexF64}(undef, ndips)

    # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
                # println(index)
                positions[index] = [xvec[kkk], yvec[kk], zvec[k]]
                magnitudes[index] = mag
                # println(pos)


            end
        end
    end
    return HertzArray(positions, [complex.([0.0, 0.0, 1.0]) for k = 1:length(positions)], magnitudes, k0)

end


using LinearAlgebra


Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = generate_AUTdips(collect(-0.25λ:λ/4:0λ), collect(-0.5λ:λ/4:0.5λ), collect(-1λ:λ/4:1λ), k0)

pwe = changerepresentation(PlaneWaveExpansion, dipoles)
length(pwe)
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

4096
``` 


```@raw html
</details>
```

```jldoctest transferexamples ; output=false
incident_pwe = AntennaFieldRepresentations.transfer(pwe, [10.0, 0.0, 0.0])
length(incident_pwe)
# output

4096
```


### `TransferMap`s
If the same transfer will be applied multiple times to the same type of `AntennaFieldRepresentation`, it is beneficial to define a corresponding `TransferMap`. A `TransferMap` represents the linear operator which takes the coefficients of an `AntennaFieldRepresentation` as input and returns the coefficients of the corresponding transferred `AntennaFieldRepresentation`.  

Analoguous to the `transfer` command, a `TransferMap` is constructed via the constructor
```julia
T = TransferMap(aut_field::AntennaFieldRepresentation, R::AbstractVector)
```

All `TransferMap`s are a subtype of the abstract type [`OperationMap`](@ref operationmaps_linmap), i.e., they behave as liniear maps.

The transpose, adjoint, and inverse operators are obtained by applying the `transpose`, `adjoint`, or `inverse` command, respectively.
