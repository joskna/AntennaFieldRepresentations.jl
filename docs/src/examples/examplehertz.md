
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
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

90-element HertzArray{Float64, ComplexF64}:
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 6.123233995736766e-17 + 1.0im
                   1.0 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
                       ⋮
    0.5403023058681398 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 6.123233995736766e-17 + 1.0im
                   1.0 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
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
julia> swe = changerepresentation(SphericalWaveExpansion, dipoles)
2046-element SphericalWaveExpansion{Radiated, SphericalCoefficients{ComplexF64}, ComplexF64}:
      372.81882940219094 - 1.138881592722361e-13im
  2.0429498849964354e-15 + 3.9968028886505635e-14im
                    -0.0 - 0.0im
       -597.139457879735 - 60.148368197526686im
      372.81882940219083 + 6.907659518598972e-14im
  -3.819306724396686e-15 - 1.4210854715202004e-14im
      16.656990543598212 - 372.327084704324im
  1.9012569296705806e-14 + 2.3322611893110605e-16im
  1.7763568394002505e-14 - 1.9884095375339e-14im
   1.183004450936826e-13 + 534.9185258328564im
                         ⋮
   5.607793188565244e-37 - 2.298091455857572e-36im
   7.931019520425135e-36 - 1.6456920911512025e-37im
  -2.388316693605139e-23 - 4.8533758573512445e-22im
   8.964487736955217e-23 - 3.526518775799417e-38im
 -1.6031609994412303e-43 - 4.448511433893094e-37im
  -3.101639572066862e-38 + 3.019091936193961e-39im
   2.457921474866289e-25 + 1.546012962146155e-24im
  -2.376296086859896e-24 - 8.364595896301584e-39im
  -4.757754865929577e-41 + 6.357635080703122e-40im
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

# output

1056-element SphericalWaveExpansion{Radiated, SphericalCoefficients{ComplexF64}, ComplexF64}:
      372.81882940219094 - 1.138881592722361e-13im
  2.0429498849964354e-15 + 3.9968028886505635e-14im
                    -0.0 - 0.0im
       -597.139457879735 - 60.148368197526686im
      372.81882940219083 + 6.907659518598972e-14im
  -3.819306724396686e-15 - 1.4210854715202004e-14im
      16.656990543598212 - 372.327084704324im
  1.9012569296705806e-14 + 2.3322611893110605e-16im
  1.7763568394002505e-14 - 1.9884095375339e-14im
   1.183004450936826e-13 + 534.9185258328564im
                         ⋮
  2.9888361977013695e-26 - 3.1082956646475002e-25im
  -7.875987277726582e-26 - 9.033242708155426e-28im
 -4.0834764403870996e-26 + 8.074125823881903e-12im
 -3.4369270330878742e-12 - 3.8497293421033457e-13im
  -5.553580772755923e-28 + 1.342658056012195e-26im
 -1.0775839965319021e-27 + 9.992108129849917e-30im
 -1.8698468464299812e-28 + 8.63390384807314e-14im
   8.004313629882451e-14 + 9.887552298066218e-15im
  -9.269115636346889e-30 - 3.886369606789213e-29im
```
```jldoctest changeexamples ; output=false
julia> equivalentorder(swe_ϵ)
22
```
or by specifying the modeorder `L` directly (in this case the estimate `ϵ` does not have any effect)
```jldoctest changeexamples ; output=false
swe_L = changerepresentation(SphericalWaveExpansion, dipoles, L=15)

# output

510-element SphericalWaveExpansion{Radiated, SphericalCoefficients{ComplexF64}, ComplexF64}:
      372.81882940219094 - 1.138881592722361e-13im
  2.0429498849964354e-15 + 3.9968028886505635e-14im
                    -0.0 - 0.0im
       -597.139457879735 - 60.148368197526686im
      372.81882940219083 + 6.907659518598972e-14im
  -3.819306724396686e-15 - 1.4210854715202004e-14im
      16.656990543598212 - 372.327084704324im
  1.9012569296705806e-14 + 2.3322611893110605e-16im
  1.7763568394002505e-14 - 1.9884095375339e-14im
   1.183004450936826e-13 + 534.9185258328564im
                         ⋮
  1.7686351878031282e-20 - 1.9956096237311316e-18im
   4.354432321711454e-19 - 7.030373462210693e-20im
    -3.19668544933544e-5 - 8.625495268927583e-5im
    -8.45476048738482e-6 - 3.3085830471794484e-20im
   2.437554500040577e-21 + 2.3505164286306834e-20im
 -2.5844742680964855e-20 + 4.208694956669805e-21im
    8.782399048595974e-7 + 3.931351355402204e-6im
  -1.1150817611979755e-6 - 4.4032156280040725e-21im
   8.056704687920506e-23 + 3.8422539453088135e-22im
```
```jldoctest changeexamples ; output=false
julia> equivalentorder(swe_L)
15
```

---

### Conversion into a `PlaneWaveExpansion{Radiated}`
```jldoctest changeexamples ; output=false
julia> pwe = changerepresentation(PlaneWaveExpansion, dipoles)
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
julia> mlfmm = changerepresentation(MLFMMSource, dipoles)
90-element MLFMMSource{HertzArray{Float64, ComplexF64}, AntennaFieldRepresentations.θϕResampleMap{AntennaFieldRepresentations.LocalθResampleMap{GaussLegendreθRegularϕSampling, GaussLegendreθRegularϕSampling, 8, Float64}, AntennaFieldRepresentations.LocalϕResampleMap{GaussLegendreθRegularϕSampling, GaussLegendreθRegularϕSampling, 8, Float64}, GaussLegendreθRegularϕSampling, GaussLegendreθRegularϕSampling, 8, 8, Float64}, GaussLegendreθRegularϕSampling, ComplexF64, AntennaFieldRepresentations.MLFMMTree{AntennaFieldRepresentations.BoxData{3, Float64}, 3, Float64}}:
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 6.123233995736766e-17 + 1.0im
                   1.0 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
                       ⋮
    0.5403023058681398 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 6.123233995736766e-17 + 1.0im
                   1.0 + 0.0im
 5.373643377032895e-17 + 0.8775825618903728im
    0.8775825618903728 + 0.0im
 3.308397447266758e-17 + 0.5403023058681398im
    0.5403023058681398 + 0.0im
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
julia> E1 = efield(dipoles, [5.0,5.0,5.0])
3-element Vector{ComplexF64}:
 -223.88429354059383 + 53.33206681034357im
 -221.49298023037227 + 57.64121750370682im
  445.30643102831215 - 120.7344653803177im
```

```jldoctest changeexamples ; output=false
julia> E2 = efield(swe, [5.0,5.0,5.0])
3-element Vector{ComplexF64}:
 -223.88429354058366 + 53.33206681036077im
 -221.49298023036235 + 57.64121750372429im
  445.30643102829134 - 120.73446538035377im
```

```jldoctest changeexamples ; output=false
julia> E3 = efield(pwe, [5.0,5.0,5.0])
3-element Vector{ComplexF64}:
 -223.88429354057865 + 53.33206681035124im
 -221.49298023035888 + 57.641217503714074im
    445.306431028284 - 120.73446538033413im
```