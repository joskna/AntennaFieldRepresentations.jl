# Interpolation of `PlaneWaveExpansion`s and `SphericalFieldSampling`s

A very useful operation on a [spherically sampled data structure](@ref spheresampling) (i.e., a [`PlaneWaveExpansion`](@ref) or a [`SphericalFieldSampling`](@ref)) is interpolation (have a look at [the theory section](@ref interpolation_theory) for background information on how the interpolation works internally). For any function ``\bm{f}(\vartheta, \varphi)`` on the sphere[^1], sampled at finitely many sampling points ``(\vartheta_i, \varphi_i)``, we want to be able to find the function value at any point ``(\vartheta, \varphi)`` which does not necessarily coincide with a point from the set of previously sampled points. 

In particular, if the function ``\bm{f}(\vartheta, \varphi)`` is sampled according to a certain [`SphereSamplingStrategy`](@ref), we want to be able to resample the function according to a different `SphereSamplingStrategy`. This resampling step is often referred to as *interpolation* in other works, but in the context of `AntennaFieldRepresentations.jl`, we are a bit more careful with our wording. In `AntennafieldRepresentations.jl` 
- the term *interpolation* means to evaluate a function ``\bm{f}(\vartheta, \varphi)`` on the sphere, sampled according to some `SphereSamplingStrategy`, at arbitrary other sampling points ``(\vartheta, \varphi)``
- the term *resampling* means to find the values ``\bm{f}(\vartheta_i, \varphi_i)`` of a function on the sphere, sampled according to a certain `SphereSamplingStrategy` *B*, from the values ``\bm{f}(\vartheta_k, \varphi_k)`` of the same function, sampled according to another `SphereSamplingStrategy` *A* - i.e., after the resampling process, the spherical function ``\bm{f}(\vartheta, \varphi)`` which was sampled according to `SphereSamplingStrategy` *A*, is now sampled according to `SphereSamplingStrategy` *B*.


[^1]: The function ``\bm{f}(\vartheta, \varphi)`` is displayed with a bold ``\bm{f}`` to indicate its vector character. All spherically sampled data structures to be interpolated in `AntennaFieldRepresentations.jl` have two components representing two independent polarizations of the electromagnetic field.

## Interpolating Spherically Sampled Data at Arbitrary Sampling Points 
We can interpolate spherically sampled data at arbitrary sampling points by using the `interpolate` method or an `InterpolateMap`. 

The basic usage of the `interpolate` method is simple (expand "Setup Code" to see how the input data is generated):

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
coeffs = zeros(ComplexF64, 30) # coefficients up to order ℓ=3
coeffs[1:16] = ComplexF64.(collect(1:16)) # populate coefficients up to order ℓ=2
sph_coefficients= SphericalCoefficients(coeffs)
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

56-element SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{ComplexF64}, ComplexF64}:
  -7.690947962703429 - 54.94328713011183im
  1.2347590807670814 - 17.5499019503607im
 -1.9610812432357951 + 24.29339569026047im
  -13.06566259783586 + 72.77747054328808im
 -52.447776992281035 - 24.280980030520702im
 -14.529244227596982 + 22.5046204404282im
  -18.28493184608163 - 35.612772659114924im
  -91.07165819770962 - 17.517272478613442im
  -57.71036020699327 + 24.66540027378141im
  57.071696514735784 + 68.5696952617538im
                     ⋮
  -66.09585421141233 - 15.910787033134408im
   6.064508461892855 - 56.40494894306287im
   65.87172180937227 + 27.58085861292713im
  -36.86962542825946 + 77.90488798849583im
  -62.17698010887305 - 40.24280874862705im
  -44.01466528296982 - 35.001116747823716im
  12.040081792091042 - 27.076420036483743im
 -14.750689156141837 - 25.30689436288148im
   18.88086132920546 - 69.67016128115739im
``` 


```@raw html
</details>
```

```jldoctest interpolateexamples ; output=false
# pwe is a PlaneWaveExpansion{Radiated}

# Arbitrary position
θ=pi/10
ϕ= pi/7.8

# interpolate pwe at new position
newEθ, newEϕ = interpolate((θ,ϕ), pwe)

# Actual values at this position
Eθ, Eϕ = -59.44801130097685 + 58.38482519439182im, 68.0278816964276 + 75.60985071712197im

# Evaluate interpolation error
abs(Eθ - newEθ) / abs(Eθ) < 0.0055 # true
abs(Eϕ - newEϕ) / abs(Eϕ) < 0.0044  # true

# output

true
```

The default interpolation method uses (approximate) local interpolation with an interpolation order of `orderθ = 12` along θ and `orderϕ = 12` along ϕ. 
For better interpolation accuracy, we can specify the keyword arguments `orderθ` and `orderϕ`

```jldoctest interpolateexamples ; output=false
# interpolate pwe at new position with better accuracy
newEθ, newEϕ = interpolate((θ,ϕ), pwe, orderθ = 18, orderϕ = 18)

# Evaluate interpolation error
abs(Eθ - newEθ) / abs(Eθ) < 0.00055 # true
abs(Eϕ - newEϕ) / abs(Eϕ) < 0.00044  # true

# output

true
```

!!! tip
    The most effective way to reduce the interpolation error is to start with oversampled data.
    Having larger oversampling in the original data is more effectively impacting the interpolation accuracy than the interpolation order.
---  

## Resampling Spherically Sampled Data According to a New `SphereSamplingStrategy`
We can resample spherically sampled data to a new `SphereSamplingStrategy` using the [`resample`](@ref) method or a `ResampleMap`. 
```jldoctest interpolateexamples ; output=false
# define the target samplingstrategy
targetsamplingstrategy= RegularθRegularϕSampling(20, 20) 

pwe2= resample(targetsamplingstrategy, pwe)

# output
440-element PlaneWaveExpansion{Radiated, RegularθRegularϕSampling, ComplexF64}:
 -109.99025406960146 + 9.483384965974622im
  -83.34184404308296 + 2.12285509029131im
  -56.94690918302272 - 2.660562491171087im
 -31.488184256653895 - 4.398675035216028im
   -5.16088313042524 - 2.921202463945016im
   24.53545058560716 + 1.6270563805741496im
    58.6865681388036 + 8.801102640571184im
   95.41569276153045 + 17.898493868901554im
   129.5475773102929 + 28.028855932563122im
   153.2492457384852 + 38.20049957468889im
                     ⋮
  -34.33171148277943 + 57.11217972670092im
 -29.926566026928807 + 29.95132952779241im
 -26.842148567855343 + 4.78300815669802im
 -29.342336093157087 - 22.1301498203163im
   -37.3959541293728 - 53.882277913601605im
 -45.207081426575435 - 90.9638794362822im
 -44.977246064396006 - 129.52052430691433im
  -29.34618795614716 - 161.4318655423414im
  4.7295500113711215 - 177.41938955071302im

```
