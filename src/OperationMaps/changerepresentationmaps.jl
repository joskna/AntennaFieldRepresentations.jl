"""
    ChangeRepresentationMap{A, B, C} <: LinearMaps.LinearMap{C}

Supertype for linear map which corresponds to a `changerepresentation` method.

# Type Parameters
- `A <: AntennaFieldRepresentation` : Type of the original representation
- `B <: AntennaFieldRepresentation` : Type of the target representation
- `C <: Complex`
"""
abstract type ChangeRepresentationMap{
    A<:AntennaFieldRepresentation,
    B<:AntennaFieldRepresentation,
    C<:Complex,
} <: LinearMaps.LinearMap{C} end

"""
    SphericalToPlaneWaveMap{S,W,C} <: ChangeRepresentationMap{S,W,C}

Linear map representing a `changerepresentation` operation from a `SphericalWaveExpansion` into a `PlaneWaveExpansion`.

# Type Parameters
- `S <: SphericalWaveExpansion` : Type of the original representation
- `W <: PlaneWaveExpansion` : Type of the target representation
- `C <: Complex`
"""
struct SphericalToPlaneWaveMap{
    S<:SphericalWaveExpansion,
    W<:PlaneWaveExpansion,
    C<:Complex,
} <: ChangeRepresentationMap{S,W,C}
    swe::S
    pwe::W
    stm::SphericalTransmitMap
end
function Base.size(stpwm::SphericalToPlaneWaveMap)
    return size(stpwm.stm)
end

function SphericalToPlaneWaveMap(
    samplingstrategy::Y,
    swe::SphericalWaveExpansion{P,H,C},
) where {
    P<:PropagationType,
    C<:Number,
    H<:AbstractSphericalCoefficients,
    Y<:SphereSamplingStrategy,
}
    Lmax = equivalentorder(swe)
    αinc = αinc_planewave(Lmax)
    fieldsampling = SphericalFieldSampling(samplingstrategy, αinc)
    pwe = PlaneWaveExpansion{P,Y,C}(
        samplingstrategy,
        fieldsampling.S21values,
        getwavenumber(swe),
    )
    stm = SphericalTransmitMap(swe, fieldsampling)
    return SphericalToPlaneWaveMap(swe, pwe, stm)
end

function SphericalToPlaneWaveMap(
    ::Type{PlaneWaveExpansion{P,Y,C}},
    swe::SphericalWaveExpansion{P,H,C},
) where {
    P<:PropagationType,
    C<:Number,
    H<:AbstractSphericalCoefficients,
    Y<:SphereSamplingStrategy,
}
    Lmax = equivalentorder(swe)
    αinc = αinc_planewave(Lmax)
    samplingstrategy = _standardsampling(S, Lmax)
    fieldsampling = SphericalFieldSampling(samplingstrategy, αinc)
    pwe = PlaneWaveExpansion{P,Y,C}(
        samplingstrategy,
        fieldsampling.S21values,
        getwavenumber(swe),
    )
    stm = SphericalTransmitMap(swe, fieldsampling)
    return SphericalToPlaneWaveMap{
        SphericalWaveExpansion{P,H,C},
        PlaneWaveExpansion{P,Y,C},
        C,
    }(
        swe,
        pwe,
        stm,
    )
end

function SphericalToPlaneWaveMap(
    ::Type{PlaneWaveExpansion},
    swe::SphericalWaveExpansion{P,H,C},
) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients}
    return SphericalToPlaneWaveMap(
        PlaneWaveExpansion{P,GaussLegendreθRegularϕSampling,C},
        swe,
    )
end

function ChangeRepresentationMap(
    ::Type{W},
    swe::SphericalWaveExpansion,
) where {W<:PlaneWaveExpansion}
    return SphericalToPlaneWaveMap(Y, swe)
end

function _linearmap(crm::SphericalToPlaneWaveMap)
    return crm.stm
end

function LinearMaps._unsafe_mul!(y, crm::ChangeRepresentationMap, x::AbstractVector)
    return LinearMaps._unsafe_mul!(y, _linearmap(crm), x)
end
