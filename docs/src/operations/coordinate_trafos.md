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
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
  133.27475990211963 - 0.28489080243388526im
   565.8482906615376 - 6.398231256339348im
   901.3284816416058 - 24.954815925545im
    741.713972906048 - 37.95216783509586im
  11.312075304375407 - 0.9195509677066037im
  -993.1008273863247 + 118.0434054654145im
 -1812.5005260873133 + 295.12568063577623im
 -2109.3918985874902 + 450.356494230454im
  -1826.429917754056 + 494.96298278464366im
 -1154.8606660853657 + 387.59518014724296im
                     ⋮
  -272.6899761452977 + 1006.231333807601im
 -281.74607593264045 + 1319.6415737761977im
 -220.95371935076062 + 1356.972842795647im
 -111.42889707612927 + 937.4750679774908im
  3.7811371002184706 - 46.325525154592im
   79.22196723618877 - 1548.2412329703982im
   92.13078373338797 - 3328.1327517310583im
    56.3835761158484 - 4989.233272438586im
  13.103463727614537 - 6097.330163053298im
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

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
  -2065.636836558667 + 398.3351821827314im
 -1546.8263526723174 + 224.5462970776981im
  -606.6235144836841 + 63.00609507931752im
   344.5805267416362 - 23.840755559205146im
   874.8217611479733 - 36.215719343276945im
   814.6349017793069 - 16.77365896141341im
     384.20122266499 - 2.6545595451635116im
   33.02700131488925 + 0.007781566122229014im
   90.63796098866803 - 0.1201697867301303im
  504.00251236073603 - 4.881003038728247im
                     ⋮
   239.9718214116018 - 219.07062228740907im
  243.23660705332668 - 258.38222023040316im
  225.84768320517733 - 279.3934697135097im
  187.02387944653907 - 270.22835119137005im
  126.47281912920127 - 214.5334757535653im
   46.98041023565534 - 94.2695547369382im
  -42.81154615469475 + 102.6782393598174im
 -127.21383309541557 + 369.6533246354198im
 -185.90786016390183 + 665.477673218149im
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

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
 -1165.4333309662043 + 249.6574114731697im
 -1233.4223089290683 + 238.60904694180823im
 -1136.5062387189096 + 202.7535155891524im
  -840.6824406923872 + 142.26856962282957im
  -338.9571833270618 + 56.26529184719954im
   343.5848325839948 - 57.9652681497453im
  1141.5695539717617 - 202.47761708612987im
  1944.3704900509672 - 373.2211102197046im
   2609.439955104748 - 554.2835753991137im
  3002.3166745968733 - 716.7376997095871im
                     ⋮
  -1356.416149328104 + 2031.356068654697im
 -1486.1590070674438 + 2508.2993980042247im
  -1552.233440184155 + 2957.236938045757im
 -1581.3084006246331 + 3405.9370716292306im
  -1603.852891367265 + 3910.1790319482525im
  -1643.068762685601 + 4535.942208238367im
  -1707.162630606429 + 5330.961423796489im
  -1786.473785869213 + 6289.059075156172im
  -1855.712163574333 + 7313.250988492382im
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

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
  -3688.734759129775 + 790.1989367265008im
  -3362.954856307884 + 650.5763243027131im
 -2685.3022711746407 + 479.05975851316316im
 -1774.2351460806449 + 300.2532583019142im
  -779.1308996346388 + 129.32950888087612im
  155.18020656282283 - 26.18100591080782im
   907.7956771720526 - 161.01358954771888im
  1394.1526394696143 - 267.60631926643566im
  1580.3206481040697 - 335.68286635400693im
  1493.5923012269818 - 356.56057232550273im
                     ⋮
  -859.6401761664597 + 1287.3888446931585im
  -742.7712340874806 + 1253.6292138258125im
 -508.67693389994395 + 969.1037967020116im
 -181.02070329533197 + 389.89484179123633im
  208.72045284691364 - 508.85859799749664im
    626.319252434787 - 1729.0508527144284im
  1034.6123798578285 - 3230.796279799118im
   1390.688634634825 - 4895.754654422732im
  1646.3759583357403 - 6488.293090745811im
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

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
  -3688.733757830636 + 790.1985666767096im
 -3362.9565329223606 + 650.5731919622826im
 -2685.3103905832904 + 479.05595666858346im
  -1774.243303996005 + 300.2509884833982im
  -779.1367523582303 + 129.3282099589832im
  155.16871213039494 - 26.183206312483424im
   907.7564869835678 - 161.0150229944847im
  1394.0978512165366 - 267.60409767123804im
  1580.3153053957544 - 335.68212317274066im
  1493.5898232301715 - 356.5372809462713im
                     ⋮
  -859.6430715013522 + 1287.391658030066im
  -742.7879909686507 + 1253.6557481299897im
  -508.7185006342291 + 969.1811730820075im
 -181.04006532605126 + 389.93870839870203im
   208.7193839778715 - 508.8556721222412im
   626.3185752037141 - 1729.0421600048899im
   1034.616578172032 - 3230.804299193661im
   1390.692580538158 - 4895.764792242206im
   1646.377006100212 - 6488.296139338653im
```

In the end, the user has to find the tradeoff between a faster or a more accurate calculation
```jldoctest rotateexamples ; output=false
using LinearAlgebra

# Evaluates to 7.537201623735068e-6 :
norm(pwe_rotχθϕ_quick - pwe_rotχθϕ) / norm(pwe_rotχθϕ)

# output

7.537201623735068e-6
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
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

4096-element PlaneWaveExpansion{Radiated, GaussLegendreθRegularϕSampling, ComplexF64}:
 -260.82365999410507 + 293.0157245851745im
  -532.7817447109824 + 697.0056758719113im
  -647.6688946792898 + 991.2719587693714im
   -519.261811688616 + 935.6284332043956im
   -99.7597635656191 + 213.48050877619963im
   574.4826299887545 - 1477.3618761741482im
  1345.2330488130806 - 4222.254551033635im
  1945.1198946389677 - 7603.226538131008im
   2092.027567606773 - 10456.380731535603im
  1645.8570546416856 - 10893.056406813364im
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
</details>
```

```jldoctest transferexamples ; output=false
incident_pwe = AntennaFieldRepresentations.transfer(pwe, [10.0, 0.0, 0.0])

# output

4096-element PlaneWaveExpansion{Incident, GaussLegendreθRegularϕSampling, ComplexF64}:
 -0.44426478475401165 + 1.0595358358054188im
    4.416249896019219 - 1.4650910235220564im
  -10.242482709670822 + 1.767763438177476im
   13.612040157130226 - 2.039889296625402im
   -4.131413583124902 + 0.7004328903529873im
   -37.06477938525286 + 7.6965649215685605im
    136.0037162092041 - 34.69368234215239im
    -314.849011927716 + 96.87967408408286im
    562.7457634256359 - 204.5112581569004im
    -778.393453947278 + 327.886483456686im
                      ⋮
                  0.0 + 0.0im
                  0.0 - 0.0im
                  0.0 + 0.0im
                  0.0 - 0.0im
                  0.0 + 0.0im
                  0.0 - 0.0im
                  0.0 + 0.0im
                  0.0 - 0.0im
                  0.0 + 0.0im
```


### `TransferMap`s
If the same transfer will be applied multiple times to the same type of `AntennaFieldRepresentation`, it is beneficial to define a corresponding `TransferMap`. A `TransferMap` represents the linear operator which takes the coefficients of an `AntennaFieldRepresentation` as input and returns the coefficients of the corresponding transferred `AntennaFieldRepresentation`.  

Analoguous to the `transfer` command, a `TransferMap` is constructed via the constructor
```julia
T = TransferMap(aut_field::AntennaFieldRepresentation, R::AbstractVector)
```

All `TransferMap`s are a subtype of the abstract type [`OperationMap`](@ref operationmaps_linmap), i.e., they behave as liniear maps.

The transpose, adjoint, and inverse operators are obtained by applying the `transpose`, `adjoint`, or `inverse` command, respectively.
