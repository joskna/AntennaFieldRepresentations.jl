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
    originalrepresentation::S
    targetrepresentation::W
    lmap::SphericalTransmitMap
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
    return SphericalToPlaneWaveMap{typeof(swe),typeof(pwe),C}(swe, pwe, stm)
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
    samplingstrategy = _standardsampling(Y, Lmax)
    return SphericalToPlaneWaveMap(samplingstrategy, swe)
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
    return SphericalToPlaneWaveMap(W, swe)
end

"""
    PlaneWaveToSphericalMap{W,S,C} <: ChangeRepresentationMap{W,S,C}

Linear map representing a `changerepresentation` operation from a `SphericalWaveExpansion` into a `PlaneWaveExpansion`.

# Type Parameters
- `W <: PlaneWaveExpansion` : Type of the original representation
- `S <: SphericalWaveExpansion` : Type of the target representation
- `C <: Complex`
"""
struct PlaneWaveToSphericalMap{
    W<:PlaneWaveExpansion,
    S<:SphericalWaveExpansion,
    C<:Complex,
} <: ChangeRepresentationMap{W,S,C}
    originalrepresentation::W
    targetrepresentation::S
    lmap::InverseSphericalTransmitMap
end
function PlaneWaveToSphericalMap(
    ::Type{SphericalWaveExpansion{P,H,C}},
    pwe::PlaneWaveExpansion{P,Y,C},
) where {
    P<:PropagationType,
    C<:Number,
    H<:AbstractSphericalCoefficients,
    Y<:SphereSamplingStrategy,
}
    Lmax = equivalentorder(pwe)
    coefficients = H(zeros(C, sℓm_to_j(2, Lmax, Lmax)))
    swe = SphericalWaveExpansion(P(), coefficients, getwavenumber(pwe))
    αinc = αinc_planewave(Lmax)
    fieldsampling = SphericalFieldSampling(pwe.samplingstrategy, αinc)
    istm = InverseSphericalTransmitMap(swe, fieldsampling)
    return PlaneWaveToSphericalMap{typeof(pwe),typeof(swe),C}(pwe, swe, istm)
end
function PlaneWaveToSphericalMap(
    ::Type{SphericalWaveExpansion},
    pwe::PlaneWaveExpansion{P,Y,C},
) where {P<:PropagationType,C<:Number,Y<:SphereSamplingStrategy}
    return PlaneWaveToSphericalMap(
        SphericalWaveExpansion{P,SphericalCoefficients{C},C},
        pwe,
    )
end
function ChangeRepresentationMap(
    ::Type{S},
    pwe::PlaneWaveExpansion,
) where {S<:SphericalWaveExpansion}
    return PlaneWaveToSphericalMap(S, pwe)
end

function _inversetype(
    ::Type{PlaneWaveToSphericalMap{W,S,C}},
) where {W<:PlaneWaveExpansion,S<:SphericalWaveExpansion,C<:Complex}
    return SphericalToPlaneWaveMap{S,W,C}
end
function _inversetype(
    ::Type{SphericalToPlaneWaveMap{S,W,C}},
) where {W<:PlaneWaveExpansion,S<:SphericalWaveExpansion,C<:Complex}
    return PlaneWaveToSphericalMap{W,S,C}
end
#################################################################################
#
#################################################################################

function _linearmap(crm::ChangeRepresentationMap)
    return crm.lmap
end
function Base.size(crm::ChangeRepresentationMap)
    return size(crm.lmap)
end
function LinearMaps._unsafe_mul!(y, crm::ChangeRepresentationMap, x::AbstractVector)
    return LinearMaps._unsafe_mul!(y, _linearmap(crm), x)
end
function LinearMaps._unsafe_mul!(
    y,
    crm_ad::LinearMaps.AdjointMap{ChangeRepresentationMap},
    x::AbstractVector,
)
    return LinearMaps._unsafe_mul!(y, adjoint(_linearmap(crm_ad.lmap)), x)
end
function LinearMaps._unsafe_mul!(
    y,
    crm_tr::LinearMaps.TransposeMap{ChangeRepresentationMap},
    x::AbstractVector,
)
    return LinearMaps._unsafe_mul!(y, transpose(_linearmap(crm_tr.lmap)), x)
end
function inverse(crm::M) where {M<:ChangeRepresentationMap}
    T = _inversetype(M)
    return T(crm.targetrepresentation, crm.originalrepresentation, inverse(crm.lmap))
end

###################################################################
# function changerepresentation(
#     T::Type{W},
#     swe::SphericalWaveExpansion{P,H,C},
# ) where {W<:PlaneWaveExpansion, P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients}
# map= ChangeRepresentationMap(
#     T,
#     swe
# )

# map.targetrepresentation .= map * deepcopy(asvector(swe))

# return map.targetrepresentation
# end
function changerepresentation(
    T::Type{A},
    originalrepresentation::A2,
) where {A<:AntennaFieldRepresentation, A2<:AntennaFieldRepresentation}
map= ChangeRepresentationMap(
    T,
    originalrepresentation
)

map.targetrepresentation .= map * deepcopy(asvector(originalrepresentation))

return map.targetrepresentation
end

# function changerepresentation(
#     ::Type{PlaneWaveExpansion},
#     swe::SphericalWaveExpansion{P,H,C},
# ) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients}

# L = equivalentorder(swe)
# samplingstrategy = GaussLegendreθRegularϕSampling(L+1, 2L+2)
# _,__, θs, ϕs = weightsandsamples(samplingstrategy)
# EθEϕ = zeros(C, length(θs), length(ϕs), 2)
# for (kθ, θ) in enumerate(θs) 
#     for (kϕ, ϕ) in enumerate(ϕs) 
#         EθEϕ[kθ, kϕ, 1], EθEϕ[kθ, kϕ, 2] = farfield(swe, (θ, ϕ))
#     end
# end
# return PlaneWaveExpansion{P,GaussLegendreθRegularϕSampling,C}(samplingstrategy, EθEϕ, getwavenumber(swe))
# end

function changerepresentation(
    ::Type{PlaneWaveExpansion},
    dipoles::DipoleArray{Pdip,E,T,C};
    ϵ = 1e-7
) where {C,Pdip,E,T}
L = equivalentorder(dipoles; ϵ = ϵ)
samplingstrategy = GaussLegendreθRegularϕSampling(L+1, 2L+2)
_,__, θs, ϕs = weightsandsamples(samplingstrategy)
EθEϕ = zeros(C, length(θs), length(ϕs), 2)
for (kθ, θ) in enumerate(θs) 
    for (kϕ, ϕ) in enumerate(ϕs) 
        EθEϕ[kθ, kϕ, 1], EθEϕ[kθ, kϕ, 2] = farfield(dipoles, (θ, ϕ))
    end
end
return PlaneWaveExpansion{Pdip,GaussLegendreθRegularϕSampling,C}(samplingstrategy, EθEϕ, getwavenumber(dipoles))

end 
