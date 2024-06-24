abstract type DipoleArray <: AntennaFieldRepresentation end

struct HertzDipoleArray{T,C} <: DipoleArray
    dipoles::Vector{HertzDipole{T,C}}
end

struct FitzgeraldDipoleArray{T,C} <: DipoleArray
    dipoles::Vector{FitzgeraldDipole{T,C}}
end
