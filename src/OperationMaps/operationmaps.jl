"""
    OperationMap{C<: Complex, A::AntennaFieldRepresentation} <: LinearMaps.LinearMap{C}

Supertype for function-like object which corresponds to a certain method.
"""
abstract type OperationMap{C<:Complex, A<:AntennaFieldRepresentation} <: LinearMaps.LinearMap{C} end
Base.length(om::OperationMap) = prod(size(om))
"""
    TransmitMap{C<: Complex, A<:AntennaFieldRepresentation, F<: FieldSampling} <: OperationMap{C, A}

Supertype for function-like object which corresponds to a `transmit` method.
"""
abstract type TransmitMap{C<:Complex, A<:AntennaFieldRepresentation, F<: FieldSampling} <: OperationMap{C, A} end

"""
    SphericalTransmitMap{Y<:SphericalWaveExpansion{Radiated}, F<:SphericalFieldSampling, C <: Complex } <: TransmitMap{C, Y, F}

Function-like object which corresponds to a `transmit` method between a `SphericalWaveExpansion{Radiated}` and a `SphericalFieldSampling`.
"""
struct SphericalTransmitMap{Y<:SphericalWaveExpansion{Radiated}, F<:SphericalFieldSampling, C <: Complex } <: TransmitMap{C, Y, F}
    swe::Y
    fs::F
    L::Integer
    Nθ::Integer 
    Lθ::Integer
    Δ::Vector{Matrix{Float64}}
    # u::Base.ReshapedArray
    # v::Base.ReshapedArray
    # v_::Base.ReshapedArray
    # v__::Base.ReshapedArray
    u::AbstractArray{C}
    v::AbstractArray{C}
    v_::AbstractArray{C}
    v__::AbstractArray{C}
    fftplanθ!::FFTW.cFFTWPlan
    fftplanϕ!::FFTW.cFFTWPlan
end
function Base.size(stm::SphericalTransmitMap)
    return length(stm.fs), length(stm.swe)
end
# function LinearMaps._unsafe_mul!(y, A::SphericalTransmitMap, x::AbstractVector)
#     # y .= asvector(A(x))
#     y .= A(x)
#     return y
# end

function SphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated, C, S},     
    fs::SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{C}, C}
    ) where {
        C<:Complex, 
        S<:AbstractSphericalCoefficients
        }
    L, Nθ, Lθ, u, v__, v_, v, S21, Δ,fftplanθ!,fftplanϕ! = _storage_fastsphrical(asvector(swe), asvector(fs.S21values), fs.samplingstrategy.Jθ, fs.samplingstrategy.Jϕ)
    return SphericalTransmitMap{SphericalWaveExpansion{Radiated, C, S}, SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{C}, C}, C}(
                    swe,
                    fs,
                    L,
                    Nθ,
                    Lθ,
                    Δ,
                    u,
                    v,
                    v_,
                    v__,
                    fftplanθ!,
                    fftplanϕ!
                    )
end


function _isfirstorder(::Type{FirstOrderSphericalCoefficients{C}}) where {C}
    return true
end
function _isfirstorder(::Type{SphericalCoefficients{C}}) where {C}
    return false
end
function (stm::SphericalTransmitMap{SphericalWaveExpansion{Radiated, C, A1}, SphericalFieldSampling{RegularθRegularϕSampling , A2, C},  C })(x::AbstractVector) where{ A1<:AbstractSphericalCoefficients, A2<:AbstractSphericalCoefficients, C<:Complex}
    firstorder= _isfirstorder(A2)
    αtoβ!(stm.swe.coefficients, x)
    return fastsphericalforward!(stm.fs.incidentcoefficients, stm.swe.coefficients, stm.fs.samplingstrategy.Jθ, stm.fs.samplingstrategy.Jϕ, stm.L, stm.Nθ, stm.Lθ, stm.u, stm.v__, stm.v_, stm.v, stm.fs.S21values, stm.Δ, stm.fftplanθ!, stm.fftplanϕ!, firstorder=firstorder)
end

function transmit(
    swe::SphericalWaveExpansion{Radiated, C, S}, 
    fs::SphericalFieldSampling{RegularθRegularϕSampling, FirstOrderSphericalCoefficients{C}, C}
    ) where {
        C<:Complex, 
        S<:AbstractSphericalCoefficients, 
        } 
        
        tmp= fastsphericalforward(
            fs.incidentcoefficients,
            αtoβ(swe.coefficients),
            fs.samplingstrategy.Jθ,
            fs.samplingstrategy.Jϕ,
            firstorder=true
        )
        fs.S21values .= tmp
        return asvector(fs.S21values)
end