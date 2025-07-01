
# MLFMMTransmitMaps
The `MLFMMTransmitMap` type calculates the transmission between a `MLFMMsource` and an `IrregularFieldSampling`.
The `MLFMMSource` is a hierarchically organized datastructure which splits the source volume into multiple box shaped domains. The source contributions of each box are aggregated and translated towards the observation locations defined in the `IrregularFieldSampling` according to the [Multilevel Fast Multipole Method](@ref mlfmm_interactions). 

The keyword arguments for the construction of a `MLFMMTransmitMap` are described in the following table

|Keyword Argument   | Short Description |
|:----------------- |:----------------- |
|`expectedaccuracy::Real` | Defines the desired approximation accuracy for all operations involving the `MLFMMsource`. Mainly affects the number of samples used to store the `PlaneWaveRepresentation`s for each box.|
|`verbose = false`| Determines if (debug) text should be printed to console for all operations involving the `MLFMMSource`.|
|`minboxlength = π / (2 * wavenumber)`| Side length of the smallest box, defaults to half a wavelength. All larger box side lengths are a ``2^n``-multiple of the `minboxlength`.|
|`orderθ=8`| Interpolation order along ``\theta`` for the `ResampleMap`s used to intrpolate from one box level to the next. Only affects `ResampleMaps` with local interpolation strategies.|
|`orderϕ=8`| Interpolation order along ``\phi`` for the `ResampleMap`s used to intrpolate from one box level to the next. Only affects `ResampleMaps` with local interpolation strategies.|
|`samplingtype::Type{SphereSamplingStrategy} = GaussLegendreθRegularϕSampling`| Type of the sampling strategy of the `PlaneWaveRepresentation`s used to represent the accumulated field contributions of each box|
|`num_bufferboxes::Integer = 1`||


MLFMMTransmitMap(
    basisfunctions, #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions)
    fieldsampling::IrregularFieldSampling,
    wavenumber::T;
    expectedaccuracy::T = T(1e-3),
    verbose = false,
    minhalfsize = π / (2 * wavenumber),
    orderθ = 8,
    orderϕ = 8,
    samplingtype::Type{S} = GaussLegendreθRegularϕSampling,
    
    transfertype::Type{AT} = PlannedTransfer{samplingtype,Complex{T},T},
    mintranslationlevel::Integer = 0,
) where {T<:Real,S<:SphereSamplingStrategy,AT<:AbstractTransfer}
    C = Complex{T}
