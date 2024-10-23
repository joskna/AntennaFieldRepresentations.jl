"""
    OperationMap{C<: Complex, A::AntennaFieldRepresentation} <: LinearMaps.LinearMap{C}

Supertype for function-like object which corresponds to a certain method.
"""
abstract type OperationMap{C<:Complex,A<:AntennaFieldRepresentation} <:
              LinearMaps.LinearMap{C} end
Base.length(om::OperationMap) = prod(size(om))
"""
    TransmitMap{C<: Complex, A<:AntennaFieldRepresentation, F<: FieldSampling} <: OperationMap{C, A}

Supertype for function-like object which corresponds to a `transmit` method.
"""
abstract type TransmitMap{C<:Complex,A<:AntennaFieldRepresentation,F<:FieldSampling} <:
              OperationMap{C,A} end

include("transmitmaps.jl")
