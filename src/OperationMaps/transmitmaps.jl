
"""
    SphericalTransmitMap{Y<:SphericalWaveExpansion{Radiated}, F<:SphericalFieldSampling, C <: Complex } <: TransmitMap{C, Y, F}

Function-like object which corresponds to a `transmit` method between a `SphericalWaveExpansion{Radiated}` and a `SphericalFieldSampling`.
"""
struct SphericalTransmitMap{
    Y<:SphericalWaveExpansion{Radiated},
    F<:SphericalFieldSampling,
    C<:Complex,
} <: TransmitMap{C,Y,F}
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
    u::Array{C}
    v::Array{C}
    v_::Array{C}
    v__::Array{C}
    fftplanθ!::FFTW.cFFTWPlan
    fftplanϕ!::FFTW.cFFTWPlan
    indicesirregθ::Vector{Vector{Int}}
    cosmθ::Vector{Vector{Float64}}
    sinmθ::Vector{Vector{ComplexF64}}
    Jθoversampled::Integer
    Jϕoversampled::Integer
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
    swe::SphericalWaveExpansion{Radiated,C,S},
    fs::SphericalFieldSampling{RegularθRegularϕSampling,S2,C},
) where {C<:Complex,S<:AbstractSphericalCoefficients,S2<:AbstractSphericalCoefficients}
    firstorder = _isfirstorder(S2)
    L, Nθ, Lθ, u, v__, v_, v, S21, Δ, fftplanθ!, fftplanϕ!, Jθoversampled, Jϕoversampled =
        _storage_fastspherical(
            asvector(swe),
            asvector(fs.S21values),
            fs.samplingstrategy.Jθ,
            fs.samplingstrategy.Jϕ,
            firstorder = firstorder,
        )
    indicesirregθ = Vector{Vector{Float64}}(undef, 0)
    cosmθ = Vector{Vector{Float64}}(undef, 0)
    sinmθ = Vector{Vector{Float64}}(undef, 0)
    return SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,S},
        SphericalFieldSampling{RegularθRegularϕSampling,S2,C},
        C,
    }(
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
        fftplanϕ!,
        indicesirregθ,
        cosmθ,
        sinmθ,
        Jθoversampled,
        Jϕoversampled,
    )
end

function SphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated,C,S},
    fs::SphericalFieldSampling{
        GaussLegendreθRegularϕSampling,
        FirstOrderSphericalCoefficients{C},
        C,
    },
) where {C<:Complex,S<:AbstractSphericalCoefficients}

    θweights, ϕweights, θs, ϕs = weightsandsamples(fs.samplingstrategy)

    L, u, v_, v, S21, cosmθ, sinmθ, indicesirregθ, Δ, fftplanϕ, Jϕoversampled =
        _storage_fastspherical_irregularθ(asvector(swe), θs, fs.samplingstrategy.Jϕ)


    return SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,S},
        SphericalFieldSampling{
            GaussLegendreθRegularϕSampling,
            FirstOrderSphericalCoefficients{C},
            C,
        },
        C,
    }(
        swe,
        fs,
        L,
        fs.samplingstrategy.Nθ,
        0,
        Δ,
        u,
        v,
        v_,
        v_,
        fftplanϕ,
        fftplanϕ,
        indicesirregθ,
        cosmθ,
        sinmθ,
        0,
        Jϕoversampled,
    )
end


function _isfirstorder(::Type{FirstOrderSphericalCoefficients{C}}) where {C}
    return true
end
function _isfirstorder(::Type{SphericalCoefficients{C}}) where {C}
    return false
end
function (
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,A1},
        SphericalFieldSampling{RegularθRegularϕSampling,A2,C},
        C,
    }
)(
    x::AbstractVector,
) where {A1<:AbstractSphericalCoefficients,A2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(A2)
    αtoβ!(stm.swe, x)

    return fastsphericalforward!(
        stm.fs.incidentcoefficients,
        stm.swe.coefficients,
        stm.fs.samplingstrategy.Jθ,
        stm.fs.samplingstrategy.Jϕ,
        stm.Jθoversampled,
        stm.Jϕoversampled,
        stm.L,
        stm.Nθ,
        stm.Lθ,
        stm.u,
        stm.v__,
        stm.v_,
        stm.v,
        stm.fs.S21values,
        stm.Δ,
        stm.fftplanθ!,
        stm.fftplanϕ!,
        firstorder = firstorder,
    )

end
function (
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,A1},
        SphericalFieldSampling{GaussLegendreθRegularϕSampling,A2,C},
        C,
    }
)(
    x::AbstractVector,
) where {A1<:AbstractSphericalCoefficients,A2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(A2)
    if !(firstorder)
        throw(
            error(
                "Irregular θ angles are only supported for first order probes at the moment.",
            ),
        )
    end
    αtoβ!(stm.swe.coefficients, x)
    return fastsphericalforward!(
        stm.fs.incidentcoefficients,
        stm.swe.coefficients,
        stm.fs.samplingstrategy.Nθ,
        stm.fs.samplingstrategy.Jϕ,
        stm.Jϕoversampled,
        stm.L,
        stm.u,
        stm.v_,
        stm.v,
        stm.fs.S21values,
        stm.cosmθ,
        stm.sinmθ,
        stm.indicesirregθ,
        stm.Δ,
        stm.fftplanϕ!,
    )

end

function transmit(
    swe::SphericalWaveExpansion{Radiated,C,S1},
    fs::SphericalFieldSampling{RegularθRegularϕSampling,S2,C},
) where {C<:Complex,S1<:AbstractSphericalCoefficients,S2<:AbstractSphericalCoefficients}

    fs.S21values .= fastsphericalforward(
        fs.incidentcoefficients,
        αtoβ(swe.coefficients),
        fs.samplingstrategy.Jθ,
        fs.samplingstrategy.Jϕ,
        firstorder = _isfirstorder(S2),
    )
    return asvector(fs.S21values)
end
function transmit(
    swe::SphericalWaveExpansion{Radiated,C,S},
    fs::SphericalFieldSampling{
        GaussLegendreθRegularϕSampling,
        FirstOrderSphericalCoefficients{C},
        C,
    },
) where {C<:Complex,S<:AbstractSphericalCoefficients}

    θweights, ϕweights, θs, ϕs = weightsandsamples(fs.samplingstrategy)
    fs.S21values .= fastsphericalforward(
        fs.incidentcoefficients,
        αtoβ(swe.coefficients),
        θs,
        fs.samplingstrategy.Jϕ,
    )
    return asvector(fs.S21values)
end


struct InverseSphericalTransmitMap{
    Y<:SphericalWaveExpansion{Radiated},
    F<:SphericalFieldSampling,
    C<:Complex,
} <: OperationMap{C,Y}
    swe::Y
    fs::F
    L::Integer
    v::Array{C}
    vview::SubArray{C}
    vview2::SubArray{C}
    ifft_planϕ::AbstractFFTs.ScaledPlan
    ifft_planθ::AbstractFFTs.ScaledPlan
    u::Matrix{C}
    P::Vector{C}
    K::Array{C}
    βaut::Vector{C}
    αaut::Vector{C}
    Amat::Matrix{C}
    uvectmp::Vector{C}
    Δ::Vector{Matrix{Float64}}
end
function Base.size(istm::InverseSphericalTransmitMap)
    return length(istm.swe), length(istm.fs)
end

function InverseSphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated,C,S},
    fs::SphericalFieldSampling{RegularθRegularϕSampling,S2,C},
) where {C<:Complex,S<:AbstractSphericalCoefficients,S2<:AbstractSphericalCoefficients}
    if !(_isfirstorder(S2))
        error("InverseSphericalTransmitMap only defined for first-order probe samplings.")
    end
    L, v, vview, vview2, ifft_planϕ, ifft_planθ, u, P, K, βaut, αaut, Amat, uvectmp, Δ =
        _storage_fastsphericalinverse(
            fs.S21values,
            fs.samplingstrategy.Jθ,
            fs.samplingstrategy.Jϕ,
        )
    return InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,S},
        SphericalFieldSampling{
            RegularθRegularϕSampling,
            FirstOrderSphericalCoefficients{C},
            C,
        },
        C,
    }(
        swe,
        fs,
        L,
        v,
        vview,
        vview2,
        ifft_planϕ,
        ifft_planθ,
        u,
        P,
        K,
        βaut,
        αaut,
        Amat,
        uvectmp,
        Δ,
    )
end
function (
    istm::InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,C,A1},
        SphericalFieldSampling{
            RegularθRegularϕSampling,
            FirstOrderSphericalCoefficients{C},
            C,
        },
        C,
    }
)(
    S12::AbstractArray,
) where {A1<:AbstractSphericalCoefficients,C<:Complex}
    istm.fs.S21values .= S12
    fastsphericalinverse!(
        istm.fs.S21values,
        istm.fs.incidentcoefficients,
        istm.fs.samplingstrategy.Jθ,
        istm.fs.samplingstrategy.Jϕ,
        istm.L,
        istm.v,
        istm.vview,
        istm.vview2,
        istm.ifft_planϕ,
        istm.ifft_planθ,
        istm.u,
        istm.P,
        istm.K,
        istm.βaut,
        istm.αaut,
        istm.Amat,
        istm.uvectmp,
        istm.Δ,
    )
    if length(istm.αaut) < length(istm.swe.coefficients)
        fill!(istm.swe.coefficients, zero(C))
        istm.swe.coefficients[1:length(istm.αaut)] .= istm.αaut
    else
        istm.swe.coefficients .= istm.αaut[1:length(istm.swe.coefficients)]
    end

    return istm.swe.coefficients
end

function inverse(stm::SphericalTransmitMap)
    return InverseSphericalTransmitMap(stm.swe, stm.fs)
end
