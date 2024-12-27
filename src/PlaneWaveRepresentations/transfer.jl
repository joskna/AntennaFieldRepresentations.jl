abstract type AbstractTransfer end
# TODO: InterpolatedTransfer

struct OnTheFlyTransfer{L,T} <: AbstractTransfer where {L<:Integer,T<:Real}
    R::SVector{3,T}
    k0::T
end
struct PlannedTransfer{C} <: AbstractTransfer where {C<:Complex}
    L::Int
    transfermatrix::Matrix{C}
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
function _initialize_plannedtransfer!(
    Pℓstorage::Vector{T},
    Rin::AbstractVector,
    k0::T,
    L::Integer,
) where {T<:Real}
    R = SVector{3,T}(Rin)
    d = cdist(R)
    kd = (k0 * d)
    h2 = collectsphericalHankel2(L + 1, kd)

    Pℓ = view(Pℓstorage, 1:L+1)
    w, θvec, ϕvec = samplingrule(L)
    transfermatrix = Matrix{Complex{T}}(undef, length(θvec), length(ϕvec))
    Rhat = SVector{3,T}(real(R) / norm(real(R)))
    sp, cp = sin.(ϕvec), cos.(ϕvec)
    st, ct = sin.(θvec), cos.(θvec)
    for k in eachindex(ϕvec)
        sinp, cosp = sp[k], cp[k]
        for kk in eachindex(θvec)
            sint, cost = st[kk], ct[kk]

            er = SVector{3,T}(sint .* cosp, sint .* sinp, cost)
            fac = Complex{T}(0.0)
            Pℓ = _collectPl!(Pℓ, L, udot(er, Rhat))
            for ℓ = 0:(L)
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ.+1] .* Pℓ[ℓ.+1]
            end
            transfermatrix[kk, k] = fac .* w[kk] * π / (2 * L + 2) / Z₀
        end
    end
    return PlannedTransfer{Complex{T}}(L, transfermatrix)
end

function _initialize_transfer!(
    _::Type{PlannedTransfer{C}},
    Pℓstorage::Vector{T},
    Rin::AbstractVector,
    k0::T,
    L::Integer,
) where {T<:Real,C}
    return _initialize_plannedtransfer!(Pℓstorage, Rin, k0, L)
end

function transfer(
    pattern::P,
    tr::OnTheFlyTransfer{T},
) where {T,P<:PlaneWaveExpansion{Radiated}}
    # L=transfer.L
    # pattern.L != L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    L = pattern.L

    patternout = converttype(reciprocaltype(P), pattern)
    R = tr.R
    d = cdist(R)
    kd = (tr.k0 * d)

    h2 = collectsphericalHankel2(L + 1, kd)

    _, θvec, ϕvec = samplingrule(pattern.L)

    for k in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[k])
        for kk in eachindex(θvec)
            sint, cost = sincos(θvec[kk])
            er = [sint * cosp; sint * sinp; cost]
            fac = Complex{T}(0, 0)
            Pℓ = collectPl(patternout.L, er ⋅ real(R) / norm(real(R)))
            for ℓ = 0:L
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ.+1] .* Pℓ[ℓ.+1]
            end
            patternout.Eθ[kk, k] *= fac
            patternout.Eϕ[kk, k] *= fac
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
    tr = OnTheFlyTransfer{pattern.L,T}(SVector{3,T}(R), getwavenumber(pattern))
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
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Radiated}}

    _muladd!(_eθ(incidentfield), _eθ(farfield), tr.transfermatrix; reset = reset)
    _muladd!(_eϕ(incidentfield), _eϕ(farfield), tr.transfermatrix; reset = reset)
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


function _adjoint_translate!(
    incidentfield::P,
    farfield::F,
    tr::PlannedTransfer{C};
    reset::Bool = true,
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Incident}}
    conj!(tr.transfermatrix)

    _muladd!(_eθ(farfield), _eθ(incidentfield), tr.transfermatrix; reset = reset)
    _muladd!(_eϕ(farfield), _eϕ(incidentfield), tr.transfermatrix; reset = reset)

    conj!(tr.transfermatrix)
    return farfield
end

function _transpose_translate!(
    incidentfield::P,
    farfield::F,
    tr::PlannedTransfer{C};
    reset::Bool = true,
) where {C,F<:PlaneWaveExpansion{Radiated},P<:PlaneWaveExpansion{Incident}}

    _muladd!(_eθ(farfield), _eθ(incidentfield), tr.transfermatrix; reset = reset)
    _muladd!(_eϕ(farfield), _eϕ(incidentfield), tr.transfermatrix; reset = reset)
    return incidentfield
end
