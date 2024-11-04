# [Transmission between an `AntennaFieldRepresentation` and a `FieldSampling`](@id transmission) 
The [`transmit`](@ref) method can be used to evaluate the probe outputs of a [`FieldSampling`](@ref fieldsampling) for a given [`AntennaFieldRepresentation`](@ref) without constructing a [`TransmitMap`](@ref) first. However, using the `transmit` method is only recommended if the operation is not called multiple times. In most cases, it is advisable to construct and use a `TransmitMap` to avoid redundant allocations. 

## [TransmitMaps](@id transmissionmap)

!!! todo
    Coming soon ...
---