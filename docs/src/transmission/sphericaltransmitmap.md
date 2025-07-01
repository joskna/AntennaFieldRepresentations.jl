
# SphericalTransmitMaps
The `SphericalTransmitMap` type calculates the transmission between a `SphericalWaveExpansion` and a `SphericalFieldSampling`. It implements the spherical Wacker algorithm which is described in more detail [in the electromagnetic theory section](@ref fastsphericalS12).

The keyword arguments for the construction of a `SphericalTransmitMap` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`firstorder::Bool` |If true, spherical coefficients of the incident probe field are treated as `FirstOrderSphericalCoefficients`. A change from `false` to `true` drops all entries which correspond to non-first-order modes. A change from `true` to `false` replaces non-first-order modes by zeros.|
|`samplingstrategy::SphereSamplingStrategy`|Replace the `SphereSamplingStrategy`. This effectively creates a new instance of `SphericalTransmitMap` from scratch.|
|`incidentcoefficients::AbstractVector`| Replace the spherical coefficients of the incident probe field by the entries of the input vector. If the `SphericalTransmitMap` is `firstorder`, all entries of the input vector which correspond to non-first-order modes are ignored.|
