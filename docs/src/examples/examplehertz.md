
## [Convert `HertzArray{Float64, ComplexF64}` into Other Representations](@id hertzarrayexample)


The following examples converts a `HertzArray{Float64, ComplexF64}` stored in the variable `dipoles` (expand "Setup Code" to see how the dipole array has been created) into various other field representations. 
```@raw html
<details closed><summary>Setup Code</summary>
```

```jldoctest changeexamples ; output=false
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
length(dipoles)
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

90

``` 


```@raw html
</details>
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/dipoles_change.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/dipoles_change.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Original field representation: An array of radiating Hertzian dipoles. 
  </figcaption>
</figure>
<br/>
```

---

### Conversion into a `SphericalWaveExpansion{Radiated}`
```jldoctest changeexamples ; output=false
julia> swe = changerepresentation(SphericalWaveExpansion, dipoles);
```


```@raw html
</details>
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/swe_change.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/swe_change.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Changed field representation: A radiating spherical wave expansion. 
  </figcaption>
</figure>
<br/>
```
The order of the resulting mode expansion (i.e., the largest considered mode order ``\ell``) is chosen according to an estimate for the resulting accuracy `ϵ` (which defaults to `1e-7`):
```jldoctest changeexamples ; output=false
julia> equivalentorder(swe)
31
```
We can tweak the resulting mode order by either specifying a new accuracy estimate `ϵ`
```jldoctest changeexamples ; output=false
swe_ϵ = changerepresentation(SphericalWaveExpansion, dipoles, ϵ =1e-2)
length(swe_ϵ)

# output

1056
```
```jldoctest changeexamples ; output=false
julia> equivalentorder(swe_ϵ)
22
```
or by specifying the modeorder `L` directly (in this case the estimate `ϵ` does not have any effect)
```jldoctest changeexamples ; output=false
swe_L = changerepresentation(SphericalWaveExpansion, dipoles, L=15)
length(swe_L)

# output

510
```
```jldoctest changeexamples ; output=false
julia> equivalentorder(swe_L)
15
```

---

### Conversion into a `PlaneWaveExpansion{Radiated}`
```jldoctest changeexamples ; output=false
julia> pwe = changerepresentation(PlaneWaveExpansion, dipoles);
```
```@raw html
</details>
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/pwe_change.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/pwe_change.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Changed field representation: A radiating plane wave expansion. 
  </figcaption>
</figure>
<br/>
```

--- 

### Conversion into a `MLFMMSource`
```jldoctest changeexamples ; output=false
julia> mlfmm = changerepresentation(MLFMMSource, dipoles);
```
```@raw html
</details>
```
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/mlfmm_change.png" width="800">
  <source media="(prefers-color-scheme: light)" srcset="../assets/mlfmm_change.png" width="800" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Changed field representation: A MLFMMSource. 
  </figcaption>
</figure>
<br/>
```

----
Of course, if we evaluate the electric (or magnetic) field, the results should be close to identical for any of the above representations (as long as the evaluation location is within the convergence region of the respective representation), as they are just different representations of the same electromagnetic field (for example consider the electric field evaluated at an arbitrary position `R=[5.0,5.0,5.0]`):

```jldoctest changeexamples ; output=false
julia> E1 = efield(dipoles, [5.0,5.0,5.0]);
```

```jldoctest changeexamples ; output=false
julia> E2 = efield(swe, [5.0,5.0,5.0]);
```

```jldoctest changeexamples ; output=false
julia> E3 = efield(pwe, [5.0,5.0,5.0]);
```