# Interpolation of `PlaneWaveExpansion`s and `SphericalFieldSampling`s

A very useful operation on a [spherically sampled data structure](@ref spheresampling) (i.e., a [`PlaneWaveExpansion`](@ref) or a [`SphericalFieldSampling`](@ref)) is interpolation (have a look at [the theory section](@ref interpolation_theory) for background information on how the interpolation works internally). For any function ``\bm{f}(\vartheta, \varphi)`` on the sphere[^1], sampled at finitely many sampling points ``(\vartheta_i, \varphi_i)``, we want to be able to find the function value at any point ``(\vartheta, \varphi)`` which does not necessarily coincide with a point from the set of previously sampled points. 

In particular, if the function ``\bm{f}(\vartheta, \varphi)`` is sampled according to a certain [`SphereSamplingStrategy`](@ref), we want to be able to resample the function according to a different `SphereSamplingStrategy`. This resampling step is often referred to as *interpolation* in other works, but in the context of `AntennaFieldRepresentations.jl`, we are a bit more careful with our wording. In `AntennafieldRepresentations.jl` 
- the term *interpolation* means to evaluate a function ``\bm{f}(\vartheta, \varphi)`` on the sphere, sampled according to some `SphereSamplingStrategy`, at arbitrary other sampling points ``(\vartheta, \varphi)``
- the term *resampling* means to find the values ``\bm{f}(\vartheta_i, \varphi_i)`` of a function on the sphere, sampled according to a certain `SphereSamplingStrategy` *B*, from the values ``\bm{f}(\vartheta_k, \varphi_k)`` of the same function, sampled according to another `SphereSamplingStrategy` *A* - i.e., after the resampling process, the spherical function ``\bm{f}(\vartheta, \varphi)`` which was sampled according to `SphereSamplingStrategy` *A*, is now sampled according to `SphereSamplingStrategy` *B*.


[^1]: The function ``\bm{f}(\vartheta, \varphi)`` is displayed with a bold ``\bm{f}`` to indicate its vector character. All spherically sampled data structures to be interpolated in `AntennaFieldRepresentations.jl` have two components representing two independent polarizations of the electromagnetic field.

## Interpolating Spherically Sampled Data at Arbitrary Sampling Points 
We can interpolate spherically sampled data at arbitrary sampling points by using the `interpolate` method or an `InterpolateMap`. 

The basic usage of the `interpolate` method is simple (expand "Setup Code") to see how the original input data is generated:

```@raw html
<details closed><summary>Setup Code</summary>
```

```jldoctest interpolateexamples ; output=false
#########################################
#                   Setup
#                     |
#                     V
#########################################
using AntennaFieldRepresentations
Z₀=376.730313669;                       # free-space wave impedance
f= 1.5e9;                               # frequency
λ = AntennaFieldRepresentations.c₀ / f; # wavelength
k0 = 2 * pi / λ;                        # wavenumber

# setup arbitrary AUT field
sph_coefficients= SphericalCoefficients(ComplexF64.(collect(1:16)))
swe= SphericalWaveExpansion(Radiated(),sph_coefficients, k0)

# Convert to PlaneWaveExpansion, by default sampled according to GaussLegendreθRegularϕSampling
pwe = changerepresentation(PlaneWaveExpansion, swe)


# Define spherical sampling strategy
_, Lmax, __ =j_to_sℓm(length(swe)); # maximum mode order
αinc=AntennaFieldRepresentations.αinc_dipole(1.0, Lmax, k0) # spherical coefficients of incident field Hertzian dipole at 1.0m distance
Jθ= 2*Lmax + 1;
Jϕ= 2*Lmax + 1;
samplingstrategy= RegularθRegularϕSampling(Jθ, Jϕ) 

# initialize empty spherical field sampling
fieldsampling = SphericalFieldSampling(samplingstrategy, αinc);
# fill fieldsampling with values
fieldsampling .= transmit(swe, fieldsampling)
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

30-element SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{ComplexF64}, ComplexF64}:
  -7.690947962703429 - 54.94328713011183im
  1.3224258014382286 - 2.603964695950758im
  -10.54097920620858 + 65.05011052267075im
 -60.343431601604024 - 4.843622334459157im
    36.6993918093313 + 10.483447973278812im
 -63.922433955021425 - 84.9789082601688im
 -29.603343764892355 + 51.94976389874797im
  115.75180538835582 + 44.2762296080817im
   50.80372394960887 - 91.33525673100621im
   42.04755897425326 + 36.950342131417585im
                     ⋮
  53.830129346425075 + 21.972387418512863im
   43.14366425476033 - 118.9949490354876im
  -90.70255751675808 - 57.17109771237522im
  44.788877769173325 - 42.61732035658832im
   59.17365937020349 + 68.38837423258568im
  -62.38029468630413 + 10.396681123456286im
  -26.14908056711137 - 48.311339908327234im
  17.055919698335316 + 13.23904623595748im
  -40.46301331669027 - 56.42703455873686im
``` 


```@raw html
</details>
```

## Resampling Spherically Sampled Data According to a New `SphereSamplingStrategy`
We can resample spherically sampled data to a new `SphereSamplingStrategy` using the `resample` method or an `ResampleMap`. 