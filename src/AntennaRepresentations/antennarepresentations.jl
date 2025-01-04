using BEAST
using StaticArrays
using LinearAlgebra


#################################################################
#
# Definition of Auxiliary types
#
#

"""
    PropagationType

Indicator for `Radiated`, `Absorbed`, or `Incident` type field expansions. 
"""
abstract type PropagationType end
struct Radiated <: PropagationType end
struct Absorbed <: PropagationType end
struct Incident <: PropagationType end

function _dualtype(P::Radiated)
    return Incident()
end
function _dualtype(P::Incident)
    return Radiated()
end

"""
    ElmagType

Indicator for `Electric` or `Magnetic` type field expansions. 
"""
abstract type ElmagType end
struct Electric <: ElmagType end
struct Magnetic <: ElmagType end

"""
    SphereSamplingStrategy

Supertype for sampling strategies on a sphere. 
"""
abstract type SphereSamplingStrategy end

"""
    RegularθRegularϕSampling <: SphereSamplingStrategy

A regular sampling strategy on the sphere has samples equiangularly distributed along the θ- and ϕ-coordinates.

ϕ ∈ {0, Δϕ, …, 2π - Δϕ} with Δϕ= 2π / Jϕ.
θ ∈ {0, Δθ, …, π - Δθ / 2} if Jθ is odd or θ ∈ {0, Δθ, …, π} if Jθ is even with Δθ = 2π / Jθ if Jθ.

# Fields: 
- `Jθ :: Integer`
- `Jϕ :: Integer`
"""
struct RegularθRegularϕSampling <: SphereSamplingStrategy
    Jθ::Integer
    Jϕ::Integer
end
function _countsamples(samplingstrategy::RegularθRegularϕSampling)
    Jθ = samplingstrategy.Jθ
    Jϕ = samplingstrategy.Jϕ
    return Jθ ÷ 2 + 1, Jϕ
end


"""
    GaussLegendreθRegularϕSampling <: SphereSamplingStrategy

Sampling strategy on the sphere with regular sampling along ϕ and Gauß-Legendre-Sampling along θ.

ϕ ∈ {0, Δϕ, …, 2π - Δϕ} with Δϕ= 2π / (2L + 2). 
θ ∈ {acos(xᵢ)} where xᵢ are the Nθ-point Gauß-Legendre quadrature points.

# Fields 
- `Nθ::Integer`
- `Jϕ::Integer`
"""
struct GaussLegendreθRegularϕSampling <: SphereSamplingStrategy
    Nθ::Integer
    Jϕ::Integer
end
function _countsamples(samplingstrategy::GaussLegendreθRegularϕSampling)
    return samplingstrategy.Nθ, samplingstrategy.Jϕ
end

"""
    _standardsampling(Type{<:SphereSamplingStrategy}, Lmax::Integer)

Return the standard `SphereSamplingStrategy` of specified type for a `PlaneWaveExpansion` with `equivalentorder=Lmax`.
"""
function _standardsampling(::Type{RegularθRegularϕSampling}, Lmax::Integer)
    Jθ = 2Lmax + 1
    Jϕ = 2Lmax + 2
    return RegularθRegularϕSampling(Jθ, Jϕ)
end
function _standardsampling(::Type{GaussLegendreθRegularϕSampling}, Lmax::Integer)
    Nθ = Lmax + 1
    Jϕ = 2Lmax + 2
    return GaussLegendreθRegularϕSampling(Nθ, Jϕ)
end

"""
    weightsandsamples(samplingstrategy::SphereSamplingStrategy) -> ( θweights::Array{Float64,1}, ϕweights::Array{Float64,1}, θs::Array{Float64,1}, ϕs::Array{Float64,1} )

Return integration weights `θweights, ϕweights` and sampling points `θs, ϕs` for the `samplingstrategy`.
"""
function weightsandsamples(samplingstrategy::RegularθRegularϕSampling)
    nθ, nϕ = _countsamples(samplingstrategy)
    dϕ = 2 * pi / (samplingstrategy.Jϕ)
    dθ = 2 * pi / (samplingstrategy.Jθ)
    ϕs = dϕ * collect(0:(nϕ-1))
    θs = dθ * collect(0:(nθ-1))
    ϕweights = fill!(Vector{Float64}(undef, nϕ), dϕ)
    θweights = sin.(θs)

    return θweights, ϕweights, θs, ϕs
end
function weightsandsamples(samplingstrategy::GaussLegendreθRegularϕSampling)
    xs::Vector{Float64}, θweights::Vector{Float64} = gausslegendre(samplingstrategy.Nθ)
    θs = acos.(-xs)
    nphi = samplingstrategy.Jϕ
    dϕ = 2 * pi / (nphi)
    ϕs = dϕ * collect(0:(nphi-1))
    ϕweights = fill!(Vector{Float64}(undef, nphi), dϕ)

    return θweights, ϕweights, θs, ϕs
end
"""
    samples(samplingstrategy::SphereSamplingStrategy) -> ( θs::Array{Float64,1}, ϕs::Array{Float64,1} )

Return sampling points `θs, ϕs` for the `samplingstrategy`.
"""
function samples(samplingstrategy::RegularθRegularϕSampling)
    nθ, nϕ = _countsamples(samplingstrategy)
    dϕ = 2 * pi / (samplingstrategy.Jϕ)
    dθ = 2 * pi / (samplingstrategy.Jθ)
    ϕs = dϕ * collect(0:(nϕ-1))
    θs = dθ * collect(0:(nθ-1))

    return θs, ϕs
end
function samples(samplingstrategy::GaussLegendreθRegularϕSampling)
    xs::Vector{Float64}, θweights::Vector{Float64} = gausslegendre(samplingstrategy.Nθ)
    θs = acos.(-xs)
    nphi = samplingstrategy.Jϕ
    dϕ = 2 * pi / (nphi)
    ϕs = dϕ * collect(0:(nphi-1))

    return θs, ϕs
end

#
# End of definition of Auxiliary types
#
#
#################################################################



#################################################################
#
# Definition of AntennaFieldRepresentation Interface
#
#

"""
    AntennaFieldRepresentation{P <: PropagationType, C <: Number} <: AbstractVector{C}

Equivalent representation of the electromagnetic fields of an antenna (in a certain region of space where the representation converges).

Each instance of `AntennaFieldRepresentation` is a discretized version of the antenna fields, i.e., a collection of coefficients from which the antenna fields can be (approximately) calculated.
The coefficents of the field representation can be accessed via the interface of `AbstractVector{C}`.
In addition to the field coefficients, each instance of `AntennaFieldRepresentation` stores additional information which allows to calculate the electromagnetic fields from the coefficient vector.
In particular, the wavenumber of the radiation frequency of the antenna field is stored.
"""
abstract type AntennaFieldRepresentation{P,C} <: AbstractVector{C} end

#
#
# Definition and docstrings of prototype functions
#
#
"""
    asvector(object::Union{AntennaFieldRepresentation,FieldSampling})

Return the coefficient vector of an `AntennaFieldRepresentation` or the measurement vector of a `FieldSampling`.
"""
function asvector end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    transmit(aut_field :: AntennaFieldRepresentation, measurement_setup :: FieldSampling{C}) -> Vector{C}

Return sampled antenna field according to the sampling defined by `measurement_setup`. 
"""
function transmit end
# must be defined for each sensible combination of instances of AntennaFieldRepresentation and FieldMeasurement


"""
    rotate(aut_field :: AntennaFieldRepresentation, χ::Number, θ::Number, ϕ::Number) -> rotated_aut_field :: AntennaFieldRepresentation

Rotate the field representation where the rotation is defined by the Euler angles `χ`, `θ`, `ϕ`.

This is equivalent to the original representation being represented in a rotated coordinate system,
where the coordinate axes of the original coordinate system must be rotated around the Euler angles
`-χ`, `-θ`, `-ϕ` to get the rotated coordinate frame.
"""
function rotate(aut_field::AntennaFieldRepresentation, χ, θ, ϕ)
    return rotate!(similar(aut_field), aut_field, χ, θ, ϕ)
end

"""
    rotate!(rotated_aut_field :: AntennaFieldRepresentation, aut_field::AntennaFieldRepresentation, χ::Number, θ::Number, ϕ::Number)

In-place rotate the field representation where the rotation is defined by the Euler angles `χ`, `θ`, `ϕ`.

This is equivalent to the `aut_field` being represented in a rotated coordinate system,
where the coordinate axes of the original coordinate system must be rotated around the Euler angles
`-χ`, `-θ`, `-ϕ` to get the rotated coordinate frame.
"""
function rotate! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    transfer(aut_field::AntennaFieldRepresentation{Radiated}, R) -> incident_field::AntennaFieldRepresentation{Incident}

Transfer radiating field representation into incident field representation in translated coordinate system with new origin at R. 
"""
function transfer(aut_field::AntennaFieldRepresentation{Radiated,C}, R) where {C}
    return transfer(similar(aut_field), aut_field, R)
end

"""
    transfer!(incident_field::AntennaFieldRepresentation{Incident}, aut_field::AntennaFieldRepresentation{Radiated}, R)

Transfer radiating field representation `aut_field` into preallocated incident field representation `incident_field` in translated coordinate system with new origin at R. 
"""
function transfer! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    spatialshift(aut_field::AntennaFieldRepresentation, R::AbstractVector) -> shifted_aut_field::AntennaFieldRepresentation
    
Spatially shift `aut_field` to new location R. 

This is equivalent to the `aut_field` being represented in a translated coordinate system,
where the coordinate axes of the original coordinate system must be translated by `-R` to get the translated coordinate frame.
"""
function spatialshift(aut_field::AntennaFieldRepresentation, R)
    return spatialshift!(similar(aut_field), aut_field, R)
end

"""
    spatialshift!(shifted_aut_field::AntennaFieldRepresentation, aut_field::AntennaFieldRepresentation, R::AbstractVector)
    
Preallocated shift of `aut_field` to new location R. 

This is equivalent to the `aut_field` being represented in a translated coordinate system,
where the coordinate axes of the original coordinate system must be translated by `-R` to get the translated coordinate frame.
"""
function spatialshift! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    efield(aut_field::AntennaFieldRepresentation, R) -> [Ex;Ey;Ez]

Return the E-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R`.

See also: [`efield!`](@ref), [`ehfield`](@ref)
"""
function efield(
    aut_field::AntennaFieldRepresentation{P,C},
    R,
) where {C<:Number,P<:PropagationType}
    storage = Vector{C}(undef, 3)
    return efield!(storage, aut_field, R)
end

"""
    efield!(storage, aut_field::AntennaFieldRepresentation, R; reset = true)

Store E-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R` in preallocated storage.
If `reset=true`, storage is overwritten. If `reset=false`, the field is added to `storage`.

See also: [`efield`](@ref), [`ehfield!`](@ref)
"""
function efield! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    hfield(aut_field::AntennaFieldRepresentation, R) -> [Hx;Hy;Hz]

Return the H-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R`.

See also: [`hfield!`](@ref), [`ehfield`](@ref)
"""
function hfield(
    aut_field::AntennaFieldRepresentation{P,C},
    R,
) where {C<:Number,P<:PropagationType}
    storage = Vector{C}(undef, 3)
    return hfield!(storage, aut_field, R)
end

"""
    hfield!(storage, aut_field::AntennaFieldRepresentation, R; reset = true)

Store H-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R` in preallocated storage.
If `reset=true`, storage is overwritten. If `reset=false`, the field is added to `storage`.

See also: [`hfield`](@ref), [`ehfield!`](@ref)
"""
function hfield! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    ehfield(aut_field::AntennaFieldRepresentation, R) -> [Ex;Ey;Ez], [Hx;Hy;Hz]

Returns the E-field and H-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R`.

For some `AntennaFieldRepresentation`s, data can be shared between the calculations of electric and magnetic fields. Therefore,
calling `ehfield` may be slightly more performant than calling `efield` and `hfield` separately in some cases.

See also: [`efield`](@ref), [`hfield`](@ref), [`ehfield!`](@ref)
"""
function ehfield(
    aut_field::AntennaFieldRepresentation{P,C},
    R,
) where {C<:Number,P<:PropagationType}
    storage_efield = Vector{C}(undef, 3)
    storage_hfield = Vector{C}(undef, 3)
    return ehfield!(storage_efield, storage_hfield, aut_field, R)
end

"""
    ehfield!(storage_efield, storage_hfield, aut_field::AntennaFieldRepresentation, R; reset=true)

Store E-field and H-field vector (in cartesian coordinates) of the field representation `aut_field` at location `R` in preallocated storages.
If `reset=true`, storages are overwritten. If `reset=false`, the fields are added to `storage_efield` and `storage_hfield`.

For some `AntennaFieldRepresentation`s, data can be shared between the calculations of electric and magnetic fields. Therefore,
calling `ehfield!` may be slightly more performant than calling `efield!` and `hfield!` separately in some cases.

See also: [`efield!`](@ref), [`hfield!`](@ref), [`ehfield`](@ref)
"""
function ehfield!(storage_efield, storage_hfield, aut_field, R; reset = true)
    return efield!(storage_efield, aut_field, R; reset = reset),
    hfield!(storage_hfield, aut_field, R; reset = reset)
end

"""
    farfield(aut_field::AntennaFieldRepresentation{Radiated}, (θ, ϕ) ) -> Eθ, Eϕ
    farfield(aut_field::AntennaFieldRepresentation{Radiated}, θ, ϕ ) -> Eθ, Eϕ   

Return Eθ,Eϕ-far-field tuple for radiating field representation into direction specified by θ and ϕ.

See also: [`efield`](@ref), [`hfield`](@ref)
"""
function farfield(
    aut_field::AntennaFieldRepresentation{Radiated,C},
    θϕ::Tuple{<:Number,<:Real},
) where {C}
    θ, ϕ = θϕ
    return farfield(aut_field, θ, ϕ)
end
function farfield(
    aut_field::AntennaFieldRepresentation{Radiated,C},
    θ::Number,
    ϕ::Real,
) where {C}
    throw(MethodError())
end

"""
    getwavenumber(aut_field::AntennaFieldRepresentation)

Returns the wavenumber k=2 π / λ which is stored in the `AntennaFieldRepresentation`.
"""
function getwavenumber(aut_field::AntennaFieldRepresentation)
    return aut_field.wavenumber
end

"""
    setwavenumber!(aut_field::AntennaFieldRepresentation, val)

Sets the wavenumber k=2 π / λ which is stored in the `AntennaFieldRepresentation` to `val`.
"""
function setwavenumber! end
# must be defined for each type which implements the interface of an AntennaFieldRepresentation

"""
    changerepresentation(Tnew::Type{<:AntennaFieldRepresentation}, aut_field::AntennaFieldRepresentation)

Return an `AntennaFieldRepresentation` of type `Tnew` which represents the same antenna fields as `aut_field`.

# Examples

```julia
# Calculate spherical wave expansion from dipole collection
shifted_hdipole=HertzArray([[0.1,0.2,0.3]], [complex.([0.0, 0.0, 1.0])],[complex(1.0)], k0)
shifted_hspherical=changerepresentation(SphericalWaveExpansion{Radiated}, shifted_hdipole)
```
"""
function changerepresentation end


"""
    equivalentorder(aut_field::AntennaFieldRepresentation; ϵ= 1e-7)

Return the estimated spherical mode order L which is needed forin a spherical mode expansion to approximate the `aut_field` to the desired accuracy `ϵ`.
"""
function equivalentorder end

"""
    inverse(opmap::OperationMap)

Return inverse linear map of `opmap`.
"""
function inverse end
#
#
# End of definition of AntennaFieldRepresentation Interface
#
#################################################################

