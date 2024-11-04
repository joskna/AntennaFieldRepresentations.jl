"""
    OperationMap{A, C} <: LinearMaps.LinearMap{C}

Supertype for function-like object which corresponds to a certain method.

# Type Parameters
- `A <: AntennaFieldRepresentation`
- `C <: Complex`
"""
abstract type OperationMap{A<:AntennaFieldRepresentation,C<:Complex} <:
              LinearMaps.LinearMap{C} end
Base.length(om::OperationMap) = prod(size(om))

"""
    TransmitMap{A, F, C} <: OperationMap{A, C}

Supertype for function-like object which corresponds to a `transmit` method.

# Type Parameters
- `A <: AntennaFieldRepresentation`
- `F <: FieldSampling`
- `C <: Complex`
"""
abstract type TransmitMap{A<:AntennaFieldRepresentation,F<:FieldSampling,C<:Complex} <:
              OperationMap{A,C} end

include("transmitmaps.jl")
