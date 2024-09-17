# [Transmission between an `AntennaFieldRepresentation` and a `FieldSampling`](@id transmission) 
The `transmit` method can be used to evaluate the probe outputs of a [`FieldSampling`](@ref fieldsampling) for a given `AntennaFieldRepresentation` without constructing a `TransmissionMap` first. However, this is only recommended if the operation is not called multiple times to avoid redundant allocations. 

## [TransmissionMaps](@id transmissionmap)

!!! note
    Coming soon ...
---