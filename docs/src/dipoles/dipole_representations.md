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

# output

81×81 Matrix{ComplexF64}:
   61.1577+109.468im  -70.2165+107.089im   …  -26.8043+198.425im
  -63.4869+112.149im  -131.709-6.13702im      -189.743+78.5838im
  -132.142+6.36302im  -65.5046-118.73im        -169.08-124.623im
  -79.7925-109.678im   66.5777-122.328im       16.1907-213.553im
   46.6733-130.774im    142.28-12.5381im       190.924-104.51im
    136.59-38.5412im   93.4168+112.515im   …   195.935+100.996im
   114.525+88.6168im  -37.4855+144.683im       30.1785+220.357im
 -0.189241+147.464im  -140.953+58.0746im      -160.982+155.067im
  -116.841+93.8382im  -135.032-76.4233im       -221.01-34.4152im
   -149.07-29.4294im  -26.8456-155.244im      -109.687-193.931im
          ⋮                                ⋱          ⋮
   94.1326-206.374im   232.896+0.913078im      206.579+10.8884im
   -110.84-198.633im   136.408-188.929im       122.936-161.009im
  -226.212-24.3534im   -77.913-219.13im       -65.6822-187.091im
  -147.366+172.706im  -224.643-56.2239im   …   -189.11-43.4156im
   62.8488+217.136im  -167.549+157.688im      -132.404+136.048im
   215.551+63.1447im    45.429+223.581im       49.4276+179.052im
   165.322-149.298im   212.975+75.0683im       176.372+43.9782im
  -43.3996-216.225im   169.599-144.995im       120.335-131.063im
  -206.647-69.3895im  -43.0956-215.873im   …  -58.8232-163.996im

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

# output

40×80 Matrix{ComplexF64}:
  1.1542e-13+942.478im   1.15055e-13+946.986im  …   1.1542e-13+942.478im
     9.60106+942.429im      -8.99759+946.934im         9.60106+942.429im
     38.3321+941.698im       1.20414+946.962im         38.3321+941.698im
     85.9184+938.553im       30.5341+946.465im         85.9184+938.553im
     151.715+930.187im        78.728+943.687im         151.715+930.187im
     234.444+912.853im       145.156+935.789im  …      234.444+912.853im
     331.879+882.112im        228.56+918.99im          331.879+882.112im
     440.516+833.193im       326.726+888.809im         440.516+833.193im
     555.309+761.51im        436.156+840.43im          555.309+761.51im
     669.543+663.307im       551.797+769.221im         669.543+663.307im
            ⋮                                   ⋱
    -555.309+761.51im       -658.528+675.046im        -555.309+761.51im
    -440.516+833.193im      -543.718+771.996im        -440.516+833.193im
    -331.879+882.112im      -428.653+842.413im        -331.879+882.112im
    -234.444+912.853im      -320.014+890.133im        -234.444+912.853im
    -151.715+930.187im      -222.798+919.805im  …     -151.715+930.187im
    -85.9184+938.553im      -140.454+936.24im         -85.9184+938.553im
    -38.3321+941.698im      -75.1564+943.903im        -38.3321+941.698im
    -9.60106+942.429im      -28.1351+946.55im         -9.60106+942.429im
 -1.1542e-13+942.478im  -1.14598e-13+946.986im     -1.1542e-13+942.478im


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