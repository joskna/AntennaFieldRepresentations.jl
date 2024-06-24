abstract type FieldSampling
end
abstract type DipoleProbes  <: FieldSampling
end

struct HertzDipoleProbes{T,C} <: DipoleProbes
    dipoles::Vector{HertzDipole{T,C}}
end

struct FitzgeraldDipoleProbes{T,C} <: DipoleProbes
    dipoles::Vector{FitzgeraldDipole{T,C}}
end

abstract type SphericalSampling  <: FieldSampling
end

struct RegularSphericalSampling{T<:Union{IncidentSphericalExpansion, FirstOrder{IncidentSphericalExpansion} }} <: SphericalSampling
    incidentexpansion::T
    Jθ::Integer 
    Jϕ::Integer 
    Jχ::Integer
end
# struct RegularSphericalSampling{ } <: SphericalSampling
#     incidentexpansion::FirstOrder{IncidentSphericalExpansion} 
#     Jθ::Integer 
#     Jϕ::Integer 
# end

# function SphericalSampling(α_inc::IncidentSphericalExpansion{C} ,Jθ, Jϕ, Jχ) where{C<:Complex}
#     return SphericalSampling(α_inc.coefficients ,Jθ, Jϕ, Jχ)
# end


struct IrregularNearFieldSampling{T<:Real, C<:Complex}
    positions::Vector{SVector{3,T}}
    θrotations::Vector{T}
    ϕrotations::Vector{T}
    χrotations::Vector{T}
    probepattern::FarfieldPattern{C}
end