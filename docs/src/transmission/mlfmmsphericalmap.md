
# MLFMMSphericalTransmitMaps
The `MLFMMSphericalTransmitMap`type calculates the transmission between an `AntennaFieldRepresentation` based on equivalent currents and a `SphericalFieldSampling`.
In its implementation, the far-field pattern (i.e., a `PlaneWaveExpansion`) is calculated from the equivalent current representation first. For maximum efficiencyc the source region is hierrchically divided into boxes arranged in an octree. The overall far-field pattern is calculated by aggregating the far-fields of sources in smaller boxes in a [MLFMM-scheme](@ref mlfmm_aggregation).
In the second step, the transmission between the calculated `PlaneWaveExpansion` and the `Spherical FieldSampling` is calculated through a `PlaneWaveSphericalTransmitMap`.

The keyword arguments for the construction of a `MLFMMSphericalTransmitMap` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`firstorder::Bool` |If true, spherical coefficients of the incident probe field are treated as `FirstOrderSphericalCoefficients`. A change from `false` to `true` drops all entries which correspond to non-first-order modes. A change from `true` to `false` replaces non-first-order modes by zeros.|
|`samplingstrategy::SphereSamplingStrategy`|Replace the `SphereSamplingStrategy`. This effectively creates a new instance of `SphericalTransmitMap` from scratch.|
|`incidentcoefficients::AbstractVector`| Replace the spherical coefficients of the incident probe field by the entries of the input vector. If the `SphericalTransmitMap` is `firstorder`, all entries of the input vector which correspond to non-first-order modes are ignored.|