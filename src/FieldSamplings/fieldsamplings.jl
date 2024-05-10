abstract type AbstractFieldSampling
end

struct HertzDipoleProbes{T,C} <: AbstractFieldSampling
    dipoles::Vector{HertzDipole{T,C}}
end

struct FitzgeraldDipoleProbes{T,C} <: AbstractFieldSampling
    dipoles::Vector{FitzgeraldDipole{T,C}}
end

struct SphericalSamplingFirstOrderProbe{C} <: AbstractFieldSampling
    coefficients::Array{C,1},  
    Jθ::Integer 
    Jϕ::Integer 
    Jχ::Integer
end
function SphericalSamplingFirstOrderProbe(α_inc::IncidentSphericalExpansion{C} ,Jθ, Jϕ, Jχ) where{C<:Complex}
    return SphericalSamplingFirstOrderProbe(α_inc.coefficients ,Jθ, Jϕ, Jχ)
end

struct SphericalSamplingAnyOrderProbe{C} <: AbstractFieldSampling
    coefficients::Array{C,1},  
    Jθ::Integer 
    Jϕ::Integer 
    Jχ::Integer
end
function SphericalSamplingAnyOrderProbe(α_inc::IncidentSphericalExpansion{C} ,Jθ, Jϕ, Jχ) where{C<:Complex}
    return SphericalSamplingAnyOrderProbe(α_inc.coefficients ,Jθ, Jϕ, Jχ)
end

struct IrregularNearFieldSampling{T<:Real, C<:Complex}
    positions::Vector{SVector{3,T}}
    θrotations::Vector{T}
    ϕrotations::Vector{T}
    χrotations::Vector{T}
    probepattern::FarfieldPattern{C}
end