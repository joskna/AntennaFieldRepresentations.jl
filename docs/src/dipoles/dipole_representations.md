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
  <source media="(prefers-color-scheme: dark)" srcset="../assets/dipoles_dark.png" width="500">
  <source media="(prefers-color-scheme: light)" srcset="../assets/dipoles_light.png" width="500" >
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    The dipoles are visualized as arrows in this example.
  </figcaption>
</figure>
<br/>
```

We might be interested in the electromagnetic field - maybe its ``E_x`` component - which is radiated by these dipoles, let's say in a plane ``z=5\lambda``.
Thus, we can define an array of observation points and evaluate the electric field at these observation points as follows:
```jldoctest dipoleexamples ; output=false
Rs= [[x , y, 5λ] for x in -10λ:λ/4:10λ , y in -10λ:λ/4:10λ] # Define observation points

E = efield.(Ref(dipoles), Rs) # Evaluate field at observation points

# output

81×81 Matrix{Vector{ComplexF64}}:
 [61.1577+109.468im, -197.913-27.1792im, -294.224+155.994im]    …  [-26.8043+198.425im, 83.6013+140.414im, -229.291+111.963im]
 [-63.4869+112.149im, -82.7781-178.009im, -293.07-160.674im]       [-189.743+78.5838im, -76.5514+132.802im, -217.555-124.047im]
 [-132.142+6.36302im, 99.99-164.503im, -32.4461-333.439im]         [-169.08-124.623im, -142.423-0.198012im, -24.6991-244.152im]
 [-79.7925-109.678im, 187.984-11.5823im, 253.695-219.302im]        [16.1907-213.553im, -66.1034-112.968im, 177.07-162.247im]
 [46.6733-130.774im, 116.757+141.955im, 326.772+74.5194im]         [190.924-104.51im, 57.3073-104.031im, 231.681+38.5753im]
 [136.59-38.5412im, -41.8016+173.939im, 138.86+304.231im]       …  [195.935+100.996im, 106.156-4.00223im, 114.833+198.947im]
 [114.525+88.6168im, -157.69+72.6061im, -151.107+296.831im]        [30.1785+220.357im, 50.9517+78.4257im, -73.5691+212.542im]
 [-0.189241+147.464im, -151.193-73.0977im, -323.036+72.5547im]     [-160.982+155.067im, -36.6205+72.3336im, -202.697+87.4858im]
 [-116.841+93.8382im, -41.8068-156.406im, -263.951-195.374im]      [-221.01-34.4152im, -69.5863+2.09718im, -200.824-83.8145im]
 [-149.07-29.4294im, 84.4934-130.536im, -31.5209-323.429im]        [-109.687-193.931im, -28.0985-53.4426im, -84.0456-198.765im]
 ⋮                                                              ⋱  ⋮
 [94.1326-206.374im, -15.2717-102.019im, -202.372+105.226im]       [206.579+10.8884im, -103.949-109.028im, -143.652+188.019im]
 [-110.84-198.633im, -85.3395-54.4505im, -15.7368+222.804im]       [122.936-161.009im, -152.047+26.1481im, 79.1442+224.247im]
 [-226.212-24.3534im, -95.6086+29.9227im, 180.651+124.686im]       [-65.6822-187.091im, -61.8064+144.813im, 232.754+53.0751im]
 [-147.366+172.706im, -36.8882+92.8956im, 195.055-93.8351im]    …  [-189.11-43.4156im, 91.8653+131.237im, 167.094-171.461im]
 [62.8488+217.136im, 50.8542+86.636im, 11.895-213.793im]           [-132.404+136.048im, 162.143-11.4585im, -62.7657-231.525im]
 [215.551+63.1447im, 100.976+11.2735im, -181.528-110.332im]        [49.4276+179.052im, 72.1085-147.903im, -232.415-60.466im]
 [165.322-149.298im, 69.0915-76.7662im, -180.477+109.842im]        [176.372+43.9782im, -95.0164-136.392im, -165.493+174.152im]
 [-43.3996-216.225im, -24.7005-102.444im, 16.9824+209.891im]       [120.335-131.063im, -166.521+19.1917im, 75.3095+228.06im]
 [-206.647-69.3895im, -101.456-36.4463im, 195.721+76.7983im]    …  [-58.8232-163.996im, -60.1649+157.684im, 236.892+38.2257im]

```

The resulting field (stored as an ordinary matrix) can then be visualized, e.g., with `Makie.jl` or `Plots.jl`
```@raw html
<figure>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../assets/field_dipoles_dark.png" width="500">
  <source media="(prefers-color-scheme: light)" srcset="../assets/field_dipoles_light.png" width="500">
  <img alt="" src="" width="200">
</picture>

  <figcaption>
    Magnitude of the Eₓ-component of the radiated field in a plane at z=5λ.
  </figcaption>
</figure>
<br/>
```



[^1]: Replacing a `Radiated` type `AntennaFieldRepresentation` by an `Absorbed` one, the electromagnetic fields of the two representations are not exactly "reversed" in time, as also the sign of the magnetic field changes. To be technically correct, both types of field representations should be considered as separate solutions of Maxwell's equations with different asymptotic boundary conditions at infinity. The fields of `DipoleArray`s of `Absorbed` type are derived from the scalar Green's function ``\mathrm{e}^{\, \mathrm{j} k r} / (4 \pi r)`` (as opposed to ``\mathrm{e}^{- \mathrm{j} k r} / (4 \pi r)`` for `Radiated` representations).

[^2]: Formally the fields of `DipoleArray`s of `Incident` type are derived from the scalar "Green's function" (more of a _pseudo_ _Green's_ _function_) ``\mathrm{sin}({\mathrm{j} k r}) / (4 \pi r)``. You can see that `Incident` fields are nothing but a superposition of `Absorbed` and `Radiated` fields.