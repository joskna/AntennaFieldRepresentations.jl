# Representing Antenna Fields with Equivalent Dipole Distributions

One of the simplest ways to represent an antenna field is by a collection of electrically short (i.e., ``\ell \ll \lambda``) dipole antennas. 
Since the [radiated fields of short dipole antennas are known analytically](@ref dipole_radiated), one can simply superimpose the effects of several spatially distributed dipole antennas to approximate the radiated fields of an antenna.

In `AntennaFieldRepresentations.jl`, collections (or arrays) of electrically short dipole arrays are stored in a struct `DipoleArray{P,E,T,C}` which is a subtype of [`AntennaFieldRepresentation{P, C}`](@ref fieldrepresentation).
The type parameters have the following meaning

| Parameter                 | Short Description                                                |
| :------------------------ | :--------------------------------------------------------------- |
| `P <: PropagationType`    | Can be `Radiated`, `Absorbed`, or `Incident`                     |
| `E <: ElmagType`          | Can be `Electric` or `Magnetic`                                  |
| `T <: Real`               | Number type used in the vector defining the positions of dipoles |
| `C <: Complex`            | Element type of the coefficient vector                           |


For extra convenience, the type aliases `HertzArray{T,C} = DipoleArray{Radiated, Electric, T, C}` and `FitzgeraldArray{T, C} = DipoleArray{Radiated, Magnetic, T, C}` are introduced. Therefore, the user will mostly interact with `HertzArray`s and `FitzgeraldArray`s while the `DipoleArray` type is hidden under the hood.

## Constructors for a `DipoleArray`
To generate a `DipoleArray`, use one of the following constructors:

```julia
DipoleArray{P, E}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
DipoleArray{P, E, T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{P <: PropagationType, E <: ElmagType, C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
HertzArray{T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
FitzgeraldArray{T, C}(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
HertzArray(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```
```julia
FitzgeraldArray(positions::Vector{V1}, orientations::Vector{V2}, dipolemoments::Vector{C}, wavenumber) where{C <: Complex>, V1<: AbstractVector, V2<: AbstractVector{C}}
```

The input arguments for the costructors are

- `positions::Vector{V1}` : A vector of 3D-position vectors. Julia must be able to convert the type `V1` into an `SVector{3}`. Must have the same length as `orientations` and `dipolemoments`.
- `orientations::Vector{V2}`: A complex valued vector of 3D-orientations. Complex values account for elliptical polarizations in general. Julia must be able to convert the type `V2` into an `SVector{3}`. Must have the same length as `positions` and `dipolemoments`.
- `dipolemoments::Vector{C}`: A vector of complex values to denote the excitation of each individual dipole. Must have the same length as `positions` and `orientations`.
- `wavenumber` : Wavenumber ``\omega = 2\pi \, f``

## Dipoles with Alternative Propagation Types
Most users will probably be familiar with radiating dipoles. They correspond to the `Radiated` propagation type. 
The `Absorbed` propagation type in some sense reverses the arrow of time[^1]. Instead of radiating power away from the dipoles towards infinity, the electromagnetic fields of an `Absorbed` propagation type bring energy from infinty towards the dipole locations. 

!!! warning
    A `DipoleArray` of `Absorbed` type must not be confused with a receiving antenna!

    The `DipoleArray` is an equivalent representation of the electromagnetic fields.
    Use the [`ProbeAntenna`](@ref probeantenna) type to indicate a receiving antenna.
---

One of the main use cases of `AntennaFieldRepresentation`s of `Absorbed` type is to represent scattered fields as a superposition of `Absorbed` and `Radiated` types. The fields of `DipoleArrays` of `Absorbed` and `Radiated` type become singular at the spots where the individual dipoles are located. 

The third type of `AntennaFieldRepresentation`, i.e., the `Incident` type[^2], does not have any singularities anywhere. It can be used to represent source-free solutions of Maxwell's equations and is well suited to represent incident fields. Thus, if field representations of `Absorbed` type are not your cup of tea, you can represent any scattered field as a superposition of `Incident` and `Radiated` types.


## Dipole Examples 

Let us first create an Array of Hertzian dipoles. 
```jldoctest dipoleexamples ; output=false
using AntennaFieldRepresentations

f = 1.5e9;  # Set frequency to 1.5 GHz
λ = AntennaFieldRepresentations.c₀ / f;  # wavelength
k0 = 2 * pi / λ;  # wavenumber

positions= [[-λ, 0, 0], [0, λ/2, λ/2], [0, -λ,0]];
orientations= [complex.([0.0,0.0,1.0]), complex.([0.0,1.0,0.0]), complex.([0.1,0.0,0.0])];
dipolemoments= [ComplexF64(1.0), ComplexF64(1.0), ComplexF64(1.0)];

dipoles = HertzArray(positions, orientations, dipolemoments, k0);

# output
3-element HertzArray{Float64, ComplexF64}:
 1.0 + 0.0im
 1.0 + 0.0im
 1.0 + 0.0im

```

The resulting set of dipoles might be visualized as follows:

```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/dipoles_light.png" width="750" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The dipoles are visualized as arrows in this example.
  </figcaption>
</figure>
<br/>
```

We might be interested in the electromagnetic field - maybe its ``E_x`` component - which is radiated by these dipoles, let's say in a plane ``z=5\lambda``.
Thus, we can define an array of observation points and evaluate the electric field at these observation points as follows (the `Ref()` command ensures that this input is treated as a constant for Julia's broadcasting operator "`.`"):
```jldoctest dipoleexamples ; output=false
Rs= [[x , y, 5λ] for x in -10λ:λ/4:10λ , y in -10λ:λ/4:10λ] # Define observation points

E = efield.(Ref(dipoles), Rs) # Evaluate E-field at observation points
Ex=[e[1] for e in E] # Extract x-component of E-field

# only show two digits after comma for output:
round.(Ex, digits = 2 )

# output

81×81 Matrix{ComplexF64}:
   61.16+109.47im   -70.22+107.09im  …  -186.01+76.12im     -26.8+198.42im
  -63.49+112.15im  -131.71-6.14im       -162.38-127.96im  -189.74+78.58im
 -132.14+6.36im      -65.5-118.73im       28.83-210.06im  -169.08-124.62im
  -79.79-109.68im    66.58-122.33im      198.74-86.63im     16.19-213.55im
   46.67-130.77im   142.28-12.54im       182.94+123.89im   190.92-104.51im
  136.59-38.54im     93.42+112.51im  …    -2.98+224.36im   195.93+101.0im
  114.52+88.62im    -37.49+144.68im     -189.06+125.68im    30.18+220.36im
   -0.19+147.46im  -140.95+58.07im      -213.66-81.78im   -160.98+155.07im
 -116.84+93.84im   -135.03-76.42im        -61.7-221.11im  -221.01-34.42im
 -149.07-29.43im    -26.85-155.24im      137.18-183.72im  -109.69-193.93im
        ⋮                            ⋱                           ⋮
   94.13-206.37im    232.9+0.91im         83.28+190.4im    206.58+10.89im
 -110.84-198.63im   136.41-188.93im      200.27+35.11im    122.94-161.01im
 -226.21-24.35im    -77.91-219.13im      134.25-146.71im   -65.68-187.09im
 -147.37+172.71im  -224.64-56.22im   …   -53.72-186.9im   -189.11-43.42im
   62.85+217.14im  -167.55+157.69im     -184.12-47.6im     -132.4+136.05im
  215.55+63.14im     45.43+223.58im      -129.4+133.61im    49.43+179.05im
  165.32-149.3im    212.97+75.07im         53.9+173.8im    176.37+43.98im
   -43.4-216.22im    169.6-145.0im       175.06+32.73im    120.34-131.06im
 -206.65-69.39im     -43.1-215.87im  …   106.33-138.2im    -58.82-164.0im

```

The resulting field (stored as an ordinary matrix) can then be visualized, e.g., with `Makie.jl` or `Plots.jl`
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/field_dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/field_dipoles_light.png" width="750" >
  <img alt="" src="" width="750">
</picture>

  <figcaption>
    Magnitude of the Eₓ-component of the radiated field in a plane at z=5λ.
  </figcaption>
</figure>
<br/>
```

Furthermore, we might be interested in the far fields radiated by the dipole collection. 
Thus, we define pairs `(θ, ϕ)` of angles for the directions in which we want to evaluate the far fields and calculate the far field via the `farfield` command 
```jldoctest dipoleexamples ; output=false
θs= LinRange(0.0, π, 40)
ϕs= LinRange(0.0, 2π, 80)

directions= [(θ, ϕ) for θ in θs, ϕ in ϕs]

farfields= farfield.(Ref(dipoles), directions)

# output

40×80 Matrix{Tuple{ComplexF64, ComplexF64}}:
 (0.0-94.2478im, 1.1542e-13+942.478im)           …  (-2.82698e-29-94.2478im, 1.1542e-13+942.478im)
 (36.7299-27.5922im, 9.60106+942.429im)             (36.7299-27.5922im, 9.60106+942.429im)
 (127.858-12.349im, 38.3321+941.698im)              (127.858-12.349im, 38.3321+941.698im)
 (225.042-76.379im, 85.9184+938.553im)              (225.042-76.379im, 85.9184+938.553im)
 (272.649-210.791im, 151.715+930.187im)             (272.649-210.791im, 151.715+930.187im)
 (231.942-374.237im, 234.444+912.853im)          …  (231.942-374.237im, 234.444+912.853im)
 (96.2882-510.728im, 331.879+882.112im)             (96.2882-510.728im, 331.879+882.112im)
 (-108.233-571.614im, 440.516+833.193im)            (-108.233-571.614im, 440.516+833.193im)
 (-334.929-531.843im, 555.309+761.51im)             (-334.929-531.843im, 555.309+761.51im)
 (-534.155-395.009im, 669.543+663.307im)            (-534.155-395.009im, 669.543+663.307im)
 ⋮                                               ⋱
 (-334.929-381.152im, -555.309+761.51im)            (-334.929-381.152im, -555.309+761.51im)
 (-108.233-412.3im, -440.516+833.193im)             (-108.233-412.3im, -440.516+833.193im)
 (96.2882-343.824im, -331.879+882.112im)            (96.2882-343.824im, -331.879+882.112im)
 (231.942-200.825im, -234.444+912.853im)            (231.942-200.825im, -234.444+912.853im)
 (272.649-31.9964im, -151.715+930.187im)         …  (272.649-31.9964im, -151.715+930.187im)
 (225.042+106.639im, -85.9184+938.553im)            (225.042+106.639im, -85.9184+938.553im)
 (127.858+173.706im, -38.3321+941.698im)            (127.858+173.706im, -38.3321+941.698im)
 (36.7299+160.292im, -9.60106+942.429im)            (36.7299+160.292im, -9.60106+942.429im)
 (8.88122e-29+94.2478im, -1.1542e-13+942.478im)     (4.278e-29+94.2478im, -1.1542e-13+942.478im)

```
The output is a pair `(Fθ, Fϕ)` for each input in the `directions` array. We can retrieve the individual components via
```jldoctest dipoleexamples ; output=false
Fθ= [ff[1] for ff in farfields]
Fϕ= [ff[2] for ff in farfields]

# only show two digits after comma for output:
round.(Fϕ, digits = 2 )

# output

40×80 Matrix{ComplexF64}:
     0.0+942.48im      0.0+946.99im  …      0.0+932.01im      0.0+942.48im
     9.6+942.43im     -9.0+946.93im       28.74+931.59im      9.6+942.43im
   38.33+941.7im       1.2+946.96im       76.35+928.98im    38.33+941.7im
   85.92+938.55im    30.53+946.47im      142.24+921.37im    85.92+938.55im
  151.71+930.19im    78.73+943.69im      225.16+905.02im   151.71+930.19im
  234.44+912.85im   145.16+935.79im  …   322.93+875.44im   234.44+912.85im
  331.88+882.11im   228.56+918.99im       432.1+827.84im   331.88+882.11im
  440.52+833.19im   326.73+888.81im      547.67+757.55im   440.52+833.19im
  555.31+761.51im   436.16+840.43im      662.95+660.74im   555.31+761.51im
  669.54+663.31im    551.8+769.22im      769.68+535.13im   669.54+663.31im
        ⋮                            ⋱
 -555.31+761.51im  -658.53+675.05im     -431.73+826.12im  -555.31+761.51im
 -440.52+833.19im  -543.72+772.0im      -322.78+874.36im  -440.52+833.19im
 -331.88+882.11im  -428.65+842.41im     -225.12+904.42im  -331.88+882.11im
 -234.44+912.85im  -320.01+890.13im     -142.24+921.1im   -234.44+912.85im
 -151.71+930.19im   -222.8+919.8im   …   -76.37+928.9im   -151.71+930.19im
  -85.92+938.55im  -140.45+936.24im      -28.75+931.6im    -85.92+938.55im
  -38.33+941.7im    -75.16+943.9im        -0.01+932.03im   -38.33+941.7im
    -9.6+942.43im   -28.14+946.55im         9.6+931.97im     -9.6+942.43im
    -0.0+942.48im     -0.0+946.99im        -0.0+932.01im     -0.0+942.48im

```

The resulting  far field can then be visualized, e.g., with `Makie.jl` or `Plots.jl`
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/farfield_dipoles_dark.png" width="750">
  <source media="(prefers-color-scheme: light)" srcset="../assets/farfield_dipoles_light.png" width="750" >
  <img alt="" src="" width="750">
</picture>

  <figcaption>
    Magnitude of  of the radiated far field in dB-scale.
  </figcaption>
</figure>
<br/>
```

As it turns out, the radiated far-field is rather omni-directional. With deeper thought, this might not be too surprising because each field component is equally well excited by the three dipoles in ``x-``, ``y-``, and ``z-`` direction.


[^1]: Replacing a `Radiated` type `AntennaFieldRepresentation` by an `Absorbed` one, the electromagnetic fields of the two representations are not exactly "reversed" in time, as also the sign of the magnetic field changes. To be technically correct, both types of field representations should be considered as separate solutions of Maxwell's equations with different asymptotic boundary conditions at infinity. The fields of `DipoleArray`s of `Absorbed` type are derived from the scalar Green's function ``\mathrm{e}^{\, \mathrm{j} k r} / (4 \pi r)`` (as opposed to ``\mathrm{e}^{- \mathrm{j} k r} / (4 \pi r)`` for `Radiated` representations).

[^2]: Formally the fields of `DipoleArray`s of `Incident` type are derived from the scalar "Green's function" (more of a _pseudo_ _Green's_ _function_) ``\mathrm{sin}({\mathrm{j} k r}) / (4 \pi r)``. You can see that `Incident` fields are nothing but a superposition of `Absorbed` and `Radiated` fields.