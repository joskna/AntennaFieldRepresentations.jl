# Interpolation of `PlaneWaveExpansion`s and `SphericalFieldSampling`s

A very useful operation on a spherically sampled data structure (i.e., a [`PlaneWaveExpansion`](@ref) or a [`SphericalFieldSampling`](@ref)) is interpolation. For any function ``\bm{f}(\vartheta, \varphi)`` on the sphere[^1], sampled at finitely many sampling points ``(\vartheta_i, \varphi_i)``, we want to be able to find the function value at any point ``(\vartheta, \varphi)`` which does not necessarily coincide with a point from the set of previously sampled points. 

In particular, if the function ``\bm{f}(\vartheta, \varphi)`` is sampled according to a certain [`SphereSamplingStrategy`](@ref), we want to be able to resample the function according to a different `SphereSamplingStrategy`. This resampling step is often referred to as *interpolation* in other works, but in the context of `AntennaFieldRepresentations.jl`, 
- the term *interpolation* means to evaluate a function ``\bm{f}(\vartheta, \varphi)`` on the sphere, sampled according to some `SphereSamplingStrategy`, at arbitrary other sampling points ``(\vartheta, \varphi)``
- the term *resampling* means to find the values ``\bm{f}(\vartheta_i, \varphi_i)`` of a function on the sphere, sampled according to a certain `SphereSamplingStrategy` *B*, from the values ``\bm{f}(\vartheta_k, \varphi_k)`` of the same function, sampled to another `SphereSamplingStrategy` *A* - i.e., after the resampling process, the spherical function ``\bm{f}(\vartheta, \varphi)`` which was sampled according to `SphereSamplingStrategy` *A*, is now sampled according to `SphereSamplingStrategy` *B*.


[^1]: The function ``\bm{f}(\vartheta, \varphi)`` is displayed with a bold ``\bm{f}`` to indicate its vector character. All spherically sampled data structures to be interpolated in `AntennaFieldRepresentations.jl` have two components representing two independent polarizations of the electromagnetic field.

## Interpolating Spherically Sampled Data at Arbitrary Sampling Points 
We can interpolate spherically sampled data at arbitrary sampling points by using the [`interpolate`](@ref) method or an [`InterpolateMap`](@ref). 

The basic usage of the `interpolate` method is simple (expand "Setup Code") to see how the original spherically sampled data is generated:

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

# Convert to PlaneWaveExpansion
pwe = changerepresentation(PlaneWaveExpansion, swe)


# Define spherical sampling strategy
_, Lmax, __ =j_to_sℓm(length(swe)); # maximum mode order
αinc=AntennaFieldRepresentations.αinc_dipole(1.0, Lmax, k0) # spherical coefficients of incident field Hertzian dipole at 1.0m distance
Jθ= 2*Lmax + 1;
Jϕ= 2*Lmax + 1;
samplingstrategy= RegularθRegularϕSampling(Jθ, Jϕ) 

# generate spherical field sampling
fieldsampling = SphericalFieldSampling(samplingstrategy, αinc);
transmit(swe, fieldsampling)
# fill fieldsampling with values
fieldsampling .= transmit(swe, fieldsampling)
###########################################
#                     ^
#                     |
#                    Setup
###########################################

# output

30-element SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{ComplexF64}, ComplexF64}:
   13.98238795699226 + 39.49725085270269im
  14.750422794923304 - 4.087573135299796im
   5.263527288822104 - 36.59295297053543im
   44.74122666311479 - 3.502849904864757im
 -36.715741326569834 + 28.07360770049009im
  0.8790221063993093 - 34.703745836516944im
  13.669210819175744 - 41.662131151398434im
  16.948092042023482 - 3.561129693815487im
  18.480191279253905 + 3.06914364570392im
   -36.2931897774764 - 22.245763190455307im
                     ⋮
 -42.602313875236995 - 9.853732645985998im
  -4.866163293003308 - 17.57655707246547im
   5.272991003993714 - 16.884288930704127im
 -26.165031009333177 + 36.57807046864393im
   5.618839566755231 + 25.692922109167313im
   6.032150250275127 + 19.5064561608415im
  26.431435394774144 + 32.46022343849675im
  35.410071426923906 - 56.70695462541424im
 -36.589220689036566 - 6.864993761385901im
``` 


```@raw html
</details>
```