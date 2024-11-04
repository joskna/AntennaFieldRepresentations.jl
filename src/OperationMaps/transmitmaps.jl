
"""
    SphericalTransmitMap{S,F,C} <: TransmitMap{S,F,C}

Linear map which corresponds to a `transmit` method between a `SphericalWaveExpansion{Radiated}` and a `SphericalFieldSampling`.

# Type Parameters
- `S <: SphericalWaveExpansion{Radiated}`
- `F <: SphericalFieldSampling`
- `C <: Complex`
"""
struct SphericalTransmitMap{
    S<:SphericalWaveExpansion{Radiated},
    F<:SphericalFieldSampling,
    C<:Complex,
} <: TransmitMap{S,F,C}
    swe::S
    fs::F
    L::Integer
    Nθ::Integer
    Lθ::Integer
    Δ::Vector{Matrix{Float64}}
    # u::Base.ReshapedArray
    # v::Base.ReshapedArray
    # v_::Base.ReshapedArray
    # v__::Base.ReshapedArray
    u::Array{C,2}
    v::Array{C,3}
    v_::Array{C,3}
    v__::Array{C,3}
    fftplanθ!::FFTW.cFFTWPlan
    fftplanϕ!::FFTW.cFFTWPlan
    cosmθ::Vector{Vector{Float64}}
    sinmθ::Vector{Vector{ComplexF64}}
    Jθoversampled::Integer
    Jϕoversampled::Integer
end
function Base.size(stm::SphericalTransmitMap)
    return length(stm.fs), length(stm.swe)
end

function TransmitMap(swe::SphericalWaveExpansion, fs::SphericalFieldSampling)
    return SphericalTransmitMap(swe, fs)
end

function SphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated,H1,C},
    fs::SphericalFieldSampling{RegularθRegularϕSampling,H,C},
) where {C<:Complex,H1<:AbstractSphericalCoefficients,H<:AbstractSphericalCoefficients}
    firstorder = _isfirstorder(H)
    L, Nθ, Lθ, u, v__, v_, v, S21, Δ, fftplanθ!, fftplanϕ!, Jθoversampled, Jϕoversampled =
        _storage_fastspherical(
            asvector(swe),
            asvector(fs.S21values),
            fs.samplingstrategy.Jθ,
            fs.samplingstrategy.Jϕ,
            firstorder = firstorder,
        )
    cosmθ = Vector{Vector{Float64}}(undef, 0)
    sinmθ = Vector{Vector{Float64}}(undef, 0)
    return SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H1,C},
        SphericalFieldSampling{RegularθRegularϕSampling,H,C},
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
        cosmθ,
        sinmθ,
        Jθoversampled,
        Jϕoversampled,
    )
end

function SphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated,H,C},
    fs::SphericalFieldSampling{
        GaussLegendreθRegularϕSampling,
        FirstOrderSphericalCoefficients{C},
        C,
    },
) where {C<:Complex,H<:AbstractSphericalCoefficients}

    θweights, ϕweights, θs, ϕs = weightsandsamples(fs.samplingstrategy)

    L, u, v_, v, S21, cosmθ, sinmθ, Δ, fftplanϕ, Jϕoversampled =
        _storage_fastspherical_irregularθ(
            fs.incidentcoefficients,
            asvector(swe),
            θs,
            fs.samplingstrategy.Jϕ,
        )


    return SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
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
        cosmθ,
        sinmθ,
        0,
        Jϕoversampled,
    )
end

"""
    InverseSphericalTransmitMap{S,F,C} <: TransmitMap{S,F,C}

Linear map which corresponds to the inverse of a `SphericalTransmitMap`.

An `InverseSphericalTransmitMap` is not necessarily available for all `SphericalTransmitMap`s, but only if efficient algorithms are known.
As a fallback, it is always possible to approximate the inverse via an iterative solver.

# Type Parameters
- `S <: SphericalWaveExpansion{Radiated}`
- `F <: SphericalFieldSampling`
- `C <: Complex`
"""
struct InverseSphericalTransmitMap{
    S<:SphericalWaveExpansion{Radiated},
    F<:SphericalFieldSampling,
    C<:Complex,
} <: OperationMap{S,C}
    swe::S
    fs::F
    L::Integer
    v::Array{C}
    w::Array{C}
    vview::SubArray{C}
    vview2::SubArray{C}
    ifft_planϕ::AbstractFFTs.ScaledPlan
    ifft_planθ::AbstractFFTs.ScaledPlan
    u::Matrix{C}
    P::Vector{C}
    K::Array{C}
    cosmθ::Vector{Vector{Float64}}
    sinmθ::Vector{Vector{C}}
    dvec::Vector{C}
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
    swe::SphericalWaveExpansion{Radiated,H,C},
    fs::SphericalFieldSampling{RegularθRegularϕSampling,H2,C},
) where {C<:Complex,H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients}
    if !(_isfirstorder(H2))
        error("InverseSphericalTransmitMap only defined for first-order probe samplings.")
    end
    L, v, vview, vview2, ifft_planϕ, ifft_planθ, u, P, K, βaut, αaut, Amat, uvectmp, Δ =
        _storage_fastsphericalinverse(
            fs.S21values,
            fs.samplingstrategy.Jθ,
            fs.samplingstrategy.Jϕ,
        )
    return InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
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
        Array{ComplexF64}(undef, 0),
        vview,
        vview2,
        ifft_planϕ,
        ifft_planθ,
        u,
        P,
        K,
        Vector{Float64}(undef, 0),
        Vector{C}(undef, 0),
        Vector{C}(undef, 0),
        βaut,
        αaut,
        Amat,
        uvectmp,
        Δ,
    )
end

function InverseSphericalTransmitMap(
    swe::SphericalWaveExpansion{Radiated,H,C},
    fs::SphericalFieldSampling{GaussLegendreθRegularϕSampling,H2,C},
) where {C<:Complex,H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients}
    if !(_isfirstorder(H2))
        error("InverseSphericalTransmitMap only defined for first-order probe samplings.")
    end
    θweights, ϕweights, θs, ϕs =
        AntennaFieldRepresentations.weightsandsamples(fs.samplingstrategy)
    w, u, v, L, ifft_planϕ!, cosmθ, sinmθ, dvec, βaut, αaut, Amat, uvectmp, Δ =
        _storage_fastsphericalinverse(fs.S21values, θs)

    return InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
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
        v,
        w,
        view(Vector{ComplexF64}(undef, 0), :),
        view(Vector{ComplexF64}(undef, 0), :),
        ifft_planϕ!,
        ifft_planϕ!,
        u,
        Vector{C}(undef, 0),
        Array{C}(undef, 0),
        cosmθ,
        sinmθ,
        dvec,
        βaut,
        αaut,
        Amat,
        uvectmp,
        Δ,
    )
end


function inverse(stm::SphericalTransmitMap)
    return InverseSphericalTransmitMap(stm.swe, stm.fs)
end

function fastsphericalforward!(
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{RegularθRegularϕSampling,H2,C},
        C,
    },
) where {H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(H2)
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

function fastsphericalforward!(
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{GaussLegendreθRegularϕSampling,H2,C},
        C,
    },
) where {H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(H2)
    if !(firstorder)
        throw(
            error(
                "Irregular θ angles are only supported for first order probes at the moment.",
            ),
        )
    end
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
        stm.Δ,
        stm.fftplanϕ!,
    )
end
function fastsphericalforward_ad!(
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{RegularθRegularϕSampling,H2,C},
        C,
    },
) where {H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(H2)
    return fastsphericalforward_ad!(
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
function fastsphericalforward_ad!(
    stm::SphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{GaussLegendreθRegularϕSampling,H2,C},
        C,
    },
) where {H<:AbstractSphericalCoefficients,H2<:AbstractSphericalCoefficients,C<:Complex}
    firstorder = _isfirstorder(H2)
    if !(firstorder)
        throw(
            error(
                "Irregular θ angles are only supported for first order probes at the moment.",
            ),
        )
    end
    return fastsphericalforward_ad!(
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
        stm.Δ,
        stm.fftplanϕ!,
    )
end

function fastsphericalinverse!(
    istm::InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{
            GaussLegendreθRegularϕSampling,
            FirstOrderSphericalCoefficients{C},
            C,
        },
        C,
    },
) where {H<:AbstractSphericalCoefficients,C<:Complex}
    θweights, ϕweights, θs, ϕs = weightsandsamples(istm.fs.samplingstrategy)
    return fastsphericalinverse!(
        istm.fs.S21values,
        istm.fs.incidentcoefficients,
        θweights,
        istm.L,
        istm.w,
        istm.u,
        istm.v,
        istm.ifft_planϕ,
        istm.cosmθ,
        istm.sinmθ,
        istm.dvec,
        istm.βaut,
        istm.αaut,
        istm.Amat,
        istm.uvectmp,
        istm.Δ,
    )
end

function fastsphericalinverse!(
    istm::InverseSphericalTransmitMap{
        SphericalWaveExpansion{Radiated,H,C},
        SphericalFieldSampling{
            RegularθRegularϕSampling,
            FirstOrderSphericalCoefficients{C},
            C,
        },
        C,
    },
) where {H<:AbstractSphericalCoefficients,C<:Complex}
    return fastsphericalinverse!(
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
end

function _isfirstorder(::Type{FirstOrderSphericalCoefficients{C}}) where {C}
    return true
end
function _isfirstorder(::Type{SphericalCoefficients{C}}) where {C}
    return false
end

function transmit(
    swe::SphericalWaveExpansion{Radiated,H,C},
    fs::SphericalFieldSampling,
) where {C<:Complex,H<:AbstractSphericalCoefficients}
    stm = SphericalTransmitMap(swe, fs)
    stm.swe .= αtoβ(stm.swe)
    y = asvector(fastsphericalforward!(stm))
    stm.swe .= βtoα(stm.swe)
    return y
end

function LinearMaps._unsafe_mul!(y, stm::SphericalTransmitMap, x::AbstractVector)
    αtoβ!(stm.swe, x)
    y .= asvector(fastsphericalforward!(stm))
    return y
end

function LinearMaps._unsafe_mul!(y, istm::InverseSphericalTransmitMap, x::AbstractVector)
    istm.fs.S21values .= reshape(x, size(istm.fs.S21values))
    fastsphericalinverse!(istm)
    if length(y) > length(istm.αaut)
        fill!(y, zero(eltype(y)))
        view(y, 1:length(istm.αaut)) .= istm.αaut
    else
        y .= view(istm.αaut, 1:length(y))
    end
    return y
end

function LinearMaps._unsafe_mul!(
    y,
    stm_ad::LinearMaps.AdjointMap{C,SphericalTransmitMap{A,B,C}},
    x::AbstractVector,
) where {A,B,C}
    stm = stm_ad.lmap
    stm.fs.S21values .= reshape(x, size(stm.fs.S21values))
    y .= βtoα!(y, fastsphericalforward_ad!(stm)) ./ 4 # division by 4 because βtoα! is 4 times the adjoint of αtoβ! .
    return y
end

function LinearMaps._unsafe_mul!(
    y,
    stm_ad::LinearMaps.TransposeMap{C,SphericalTransmitMap{A,B,C}},
    x::AbstractVector,
) where {A,B,C}
    stm = stm_ad.lmap
    stm.fs.S21values .= reshape(conj.(x), size(stm.fs.S21values))
    y .= conj.(βtoα!(y, fastsphericalforward_ad!(stm)) ./ 4) # division by 4 because βtoα! is 4 times the adjoint of αtoβ! .
    return y
end
