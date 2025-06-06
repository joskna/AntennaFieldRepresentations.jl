abstract type AbstractTransfer end
# TODO: InterpolatedTransfer

struct OnTheFlyTransfer{S,T} <: AbstractTransfer where {S<:SphereSamplingStrategy,T<:Real}
    R::SVector{3,T}
    k0::T
    L::Integer
    sampling::S
end
struct PlannedTransfer{S,C,T} <:
       AbstractTransfer where {S<:SphereSamplingStrategy,C<:Complex,T<:Real}
    R::SVector{3,T}
    k0::T
    L::Integer
    sampling::S
    transfermatrix::Matrix{C}
end

function gettransfervector(transfer::AbstractTransfer)
    return transfer.R
end
function getwavenumber(transfer::AbstractTransfer)
    return transfer.k0
end
function equivalentorder(transfer::AbstractTransfer)
    return transfer.L
end
function getsampling(transfer::AbstractTransfer)
    return transfer.sampling
end
# function _initialize_plannedtransfer(Rin::AbstractVector, k0::T, L::Integer) where{T<:Real}
#     R=SVector{3,T}(Rin)
#     d = cdist(R)
#     kd = (k0 * d)
#     h2 = collectsphericalHankel2(L + 1, kd)

#     w, θvec, ϕvec = samplingrule(L)
#     transfermatrix=Array{Complex{T}}(undef, length(θvec), length(ϕvec))
#     for k in eachindex(ϕvec)
#         sinp, cosp = sin(ϕvec[k]), cos(ϕvec[k])
#         for kk in eachindex(θvec)
#             sint, cost = sin(θvec[kk]), cos(θvec[kk])
#             er = [sint .* cosp; sint .* sinp; cost]
#             fac = Complex{T}(0.0)
#             Pℓ = collectPl(L, er ⋅ real(R) / norm(real(R)))
#             for ℓ in 0:(L)
#                 fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ .+ 1] .* Pℓ[ℓ .+ 1]
#             end
#             transfermatrix[kk, k] = fac .* w[kk] * π / (2 * L + 2) / Z₀ 
#         end
#     end
#     return PlannedTransfer{L, Complex{T}}(transfermatrix)
# end
function _initialize_transfermatrix!(
    Pℓstorage::Vector{T},
    Rin::AbstractVector,
    k0::T,
    sampling::S,
    L::Integer;
    multiplyweights::Bool = false,
) where {T<:Real,S<:SphereSamplingStrategy}

    # R = SVector{3,T}(Rin)
    R = Rin
    d = cdist(R)
    kd = (k0 * d)
    h2 = collectsphericalHankel2(L + 1, kd)

    Pℓ = view(Pℓstorage, 1:L+1)
    θweights, ϕweights, θs, ϕs = weightsandsamples(sampling)

    C = Complex{T}

    transfermatrix = Matrix{C}(undef, length(θs), length(ϕs))
    # Rhat = SVector{3,T}(real(R) / norm(real(R)))
    Rhat = real.(R) / norm(real.(R))
    sp, cp = sin.(ϕs), cos.(ϕs)
    st, ct = sin.(θs), cos.(θs)
    er = zeros(T, 3)
    for k in eachindex(ϕs)
        sinp, cosp = sp[k], cp[k]
        for kk in eachindex(θs)
            sint, cost = st[kk], ct[kk]

            er .= sint .* cosp, sint .* sinp, cost
            fac = C(0.0)
            Pℓ .= _collectPl!(Pℓ, L, udot(er, Rhat))
            for ℓ = 0:(L)
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ.+1] .* Pℓ[ℓ.+1]
            end
            # transfermatrix[kk, k] = C(0, -0.5) * fac / pi * k0
            transfermatrix[kk, k] = 0.5 * fac / Z₀
            if multiplyweights
                transfermatrix[kk, k] *= θweights[kk] * ϕweights[k]
            end
            #.* w[kk] * π / (2 * L + 2) / Z₀
        end
    end
    return transfermatrix
end

function _initialize_plannedtransfer!(
    Pℓstorage::Vector{T},
    Rin::AbstractVector,
    k0::T,
    sampling::S,
    L::Integer;
    multiplyweights::Bool = false,
) where {T<:Real,S<:SphereSamplingStrategy}
    transfermatrix = _initialize_transfermatrix!(
        Pℓstorage,
        Rin,
        k0,
        sampling,
        L;
        multiplyweights = multiplyweights,
    )

    return PlannedTransfer{S,Complex{T},T}(Rin, k0, L, sampling, transfermatrix)
end

# function _initialize_transfer!(
#     _::Type{PlannedTransfer{C}},
#     Pℓstorage::Vector{T},
#     Rin::AbstractVector,
#     k0::T,
#     L::Integer,
# ) where {T<:Real,C}
#     return _initialize_plannedtransfer!(Pℓstorage, Rin, k0, L)
# end

function transfer(
    pattern::PlaneWaveExpansion{Radiated,Y,C},
    tr::OnTheFlyTransfer{T},
) where {T,Y,C}
    # L=transfer.L
    # pattern.L != L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    L = tr.L

    a, b = size(_eθ(pattern))

    patternout = PlaneWaveExpansion(
        Incident(),
        pattern.samplingstrategy,
        Matrix(_eθ(pattern)),
        Matrix(_eϕ(pattern)),
        getwavenumber(pattern),
    )

    # println(_eθ(patternout))

    R = tr.R
    d = cdist(R)
    kd = (tr.k0 * d)

    h2 = collectsphericalHankel2(L + 1, kd)


    θvec, ϕvec = samples(pattern.samplingstrategy)

    for k in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[k])
        for kk in eachindex(θvec)
            sint, cost = sincos(θvec[kk])
            er = [sint * cosp; sint * sinp; cost]
            fac = C(0)
            Pℓ = collectPl(tr.L, er ⋅ real(R) / norm(real(R)))
            for ℓ = 0:L
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ.+1] .* Pℓ[ℓ.+1]
            end
            fac = C(0, -0.25) * fac / pi * getwavenumber(pattern)
            _eθ(patternout)[kk, k] = _eθ(pattern)[kk, k] * fac
            _eϕ(patternout)[kk, k] = _eϕ(pattern)[kk, k] * fac
        end
    end
    return patternout
end


function _imaginarypowerofℓ(ℓ::I) where {I<:Integer}
    ℓmod4 = mod(ℓ, 4)
    if ℓmod4 == 0
        return Complex{float(I)}(1, 0)
    elseif ℓmod4 == 1
        return Complex{float(I)}(0, -1)
    elseif ℓmod4 == 2
        return Complex{float(I)}(-1, 0)
    elseif ℓmod4 == 3
        return Complex{float(I)}(0, 1)
    else
        return Complex{float(I)}(0, 0)
    end
end

function transfer(
    pattern::PlaneWaveExpansion{Radiated},
    R::AbstractVector{T},
) where {T<:Real}
    L = equivalentorder(pattern)
    tr = OnTheFlyTransfer{typeof(pattern.samplingstrategy),T}(
        SVector{3,T}(R),
        getwavenumber(pattern),
        L,
        pattern.samplingstrategy,
    )
    return transfer(pattern, tr)
end

function transfer(
    pattern::PlaneWaveExpansion{Radiated,Y,C},
    tr::PlannedTransfer{C},
) where {Y,C}
    return PlaneWaveExpansion{Incident,Y,C}(
        pattern.samplingstrategy,
        _eθ(pattern) .* tr.transfermatrix,
        _eϕ(pattern) .* tr.transfermatrix,
        getwavenumber(pattern),
    )
end

function transfer!(
    incidentfield::P,
    farfield::F,
    tr::PlannedTransfer{C};
    reset::Bool = true,
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Incident}}

    _muladd_or_mulreset!(
        _eθ(incidentfield),
        _eθ(farfield),
        tr.transfermatrix;
        reset = reset,
    )
    _muladd_or_mulreset!(
        _eϕ(incidentfield),
        _eϕ(farfield),
        tr.transfermatrix;
        reset = reset,
    )
    return incidentfield
end
# function translate!(incidentfield::P, farfield::F, transfer::OnTheFlyTransfer{L,C}; reset::Bool=true) where{L, C, F<:FarfieldPattern, P<:PlaneWaveSpectrum}  
#     transferplan= PlannedTransfer(transfer.R, transfer.k0, transfer.L)
#     return translate!(incidentfield, farfield, transferplan, reset=reset)
# end
function transfer!(
    incidentfield::PlaneWaveExpansion{Incident},
    farfield::PlaneWaveExpansion{Radiated},
    R::AbstractVector{T};
    reset::Bool = true,
) where {T<:Real}
    transfer = OnTheFlyTransfer{pattern.L,T}(SVector{3,T}(R), getwavenumber(farfield))
    return transfer!(incidentfield, farfield, transfer, reset = reset)
end


function _adjoint_transfer!(
    incidentfield::P,
    farfield::F,
    tr::PlannedTransfer{C};
    reset::Bool = true,
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Incident}}
    conj!(tr.transfermatrix)

    _muladd_or_mulreset!(
        _eθ(farfield),
        _eθ(incidentfield),
        tr.transfermatrix;
        reset = reset,
    )
    _muladd_or_mulreset!(
        _eϕ(farfield),
        _eϕ(incidentfield),
        tr.transfermatrix;
        reset = reset,
    )

    conj!(tr.transfermatrix)
    return farfield
end

function _transpose_transfer!(
    incidentfield::P,
    farfield::F,
    tr::PlannedTransfer{C};
    reset::Bool = true,
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Incident}}

    _muladd_or_mulreset!(
        _eθ(farfield),
        _eθ(incidentfield),
        tr.transfermatrix;
        reset = reset,
    )
    _muladd_or_mulreset!(
        _eϕ(farfield),
        _eϕ(incidentfield),
        tr.transfermatrix;
        reset = reset,
    )
    return incidentfield
end



"""
    collectPl(Lmax,x)

Return Legendre polynomials up to Lmax
"""
function collectPl(Lmax::I, x::T) where {I<:Integer,T<:Number}

    Pℓ = zeros(T, Lmax + 1)
    Pℓ[1] = one(T)
    if Lmax > 0
        Pℓ[2] = x
    end

    # use two-term recurrence relation for Pℓ in direction of increasing ℓ
    for ℓ = 2:Lmax
        Pℓ[ℓ+1] = ((2 * ℓ - 1) * x * Pℓ[ℓ] - (ℓ - 1) * Pℓ[ℓ-1]) / ℓ
    end

    return Pℓ
end
"""
    _collectPl!(Pℓstorage, Lmax,x)

Return Legendre polynomials up to Lmax with preallocated storage
"""
function _collectPl!(
    Pℓstorage::AbstractVector{T},
    Lmax::I,
    x::T,
) where {I<:Integer,T<:Number}

    # Pℓ = zeros(T, Lmax + 1)
    Pℓstorage[1] = one(T)
    if Lmax > 0
        Pℓstorage[2] = x
    end

    # use two-term recurrence relation for Pℓ in direction of increasing ℓ
    for ℓ = 2:Lmax
        # floatℓ = T(ℓ)
        floatℓ = ℓ
        Pℓstorage[ℓ+1] =
            (
                (T(2) .* floatℓ .- T(1)) .* x .* Pℓstorage[ℓ] -
                (floatℓ .- T(1)) .* Pℓstorage[ℓ-1]
            ) ./ floatℓ
    end

    return Pℓstorage
end

function collectsphericalHankel2(Lmax::Integer, kA::N) where {N<:Number}
    zℓ = zeros(complex(N), maximum((Lmax + 1, 2)))
    expfac = cis(-kA)
    zℓ[1] = 1im * expfac / kA

    if Lmax > 0
        zℓ[2] = (1im - kA) * expfac / (kA^2)
        for ℓ = 2:Lmax
            zℓ[ℓ+1] = (2 * ℓ - 1) / kA * zℓ[ℓ] - zℓ[ℓ-1]
        end
    end

    return zℓ
end
