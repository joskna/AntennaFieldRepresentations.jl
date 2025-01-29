"""
    OperationMap{C} <: LinearMaps.LinearMap{C}

Supertype for function-like object which corresponds to a certain method acting on an `AntennaFieldRepresentation`.

# Type Parameters
- `C <: Complex`
"""
abstract type OperationMap{C<:Complex} <:
              LinearMaps.LinearMap{C} end

