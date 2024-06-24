abstract type AbstractTransfer end
# TODO: InterpolatedTransfer

struct OnTheFlyTransfer{L, T} <: AbstractTransfer where{L<:Integer, T<:Real}
    R::SVector{3,T}
    k0::T
end
struct PlannedTransfer{C} <: AbstractTransfer where{C<:Complex}
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
function _initialize_plannedtransfer!(Pℓstorage::Vector{T}, Rin::AbstractVector, k0::T, L::Integer) where{T<:Real}
    R=SVector{3,T}(Rin)
    d = cdist(R)
    kd = (k0 * d)
    h2 = collectsphericalHankel2(L + 1, kd)

    Pℓ=view(Pℓstorage, 1:L+1)
    w, θvec, ϕvec = samplingrule(L)
    transfermatrix=Matrix{Complex{T}}(undef, length(θvec), length(ϕvec))
    Rhat= SVector{3,T}(real(R) / norm(real(R)))
    sp, cp = sin.(ϕvec), cos.(ϕvec)
    st, ct = sin.(θvec), cos.(θvec)
    for k in eachindex(ϕvec)
        sinp, cosp = sp[k] , cp[k]
        for kk in eachindex(θvec)
            sint, cost = st[kk] , ct[kk]
            
            er = SVector{3,T}(sint .* cosp, sint .* sinp, cost)
            fac = Complex{T}(0.0)
            Pℓ = _collectPl!(Pℓ, L, udot(er , Rhat))
            for ℓ in 0:(L)
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ .+ 1] .* Pℓ[ℓ .+ 1]
            end
            transfermatrix[kk, k] = fac .* w[kk] * π / (2 * L + 2) / Z₀ 
        end
    end
    return PlannedTransfer{Complex{T}}(L, transfermatrix)
end

function _initialize_transfer!(_::Type{PlannedTransfer{C}}, Pℓstorage::Vector{T}, Rin::AbstractVector, k0::T, L::Integer) where{T<:Real, C}
    return _initialize_plannedtransfer!(Pℓstorage, Rin, k0, L)
end

function translate(pattern::P, transfer::OnTheFlyTransfer{T}) where{T,P<:FarfieldPattern}
    # L=transfer.L
    # pattern.L != L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")

    patternout = converttype(reciprocaltype(P), pattern)
    R = transfer.R
    d = cdist(R)
    kd = (transfer.k0 * d)

    h2 = collectsphericalHankel2(L + 1, kd)

    _, θvec, ϕvec = samplingrule(pattern.L)

    for k in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[k])
        for kk in eachindex(θvec)
            sint, cost = sincos(θvec[kk])
            er = [sint * cosp; sint * sinp; cost]
            fac = Complex{T}(0, 0)
            Pℓ = collectPl(patternout.L, er ⋅ real(R) / norm(real(R)))
            for ℓ in 0:L
                fac += _imaginarypowerofℓ(ℓ) .* (2 .* ℓ .+ 1) .* h2[ℓ .+ 1] .* Pℓ[ℓ .+ 1]
            end
            patternout.Eθ[kk, k] *= fac
            patternout.Eϕ[kk, k] *= fac
        end
    end
    return patternout
end


function _imaginarypowerofℓ(ℓ::I) where{I<:Integer}
    ℓmod4 = mod(ℓ,4)
    if  ℓmod4 == 0
        return Complex{float(I)}(1,0)
    elseif ℓmod4 == 1
        return Complex{float(I)}(0,-1)
    elseif ℓmod4 == 2
        return Complex{float(I)}(-1,0)
    elseif ℓmod4 == 3
        return Complex{float(I)}(0,1)
    else return Complex{float(I)}(0,0)
    end
end

function translate(pattern::FarfieldPattern, R::AbstractVector{T}, k₀::T) where{T<:Real}
     transfer=OnTheFlyTransfer{pattern.L, T}(SVector{3,T}(R), k₀)
    return translate(pattern, transfer)
end

function translate(pattern::P, transfer::PlannedTransfer{C}) where{C, P<:FarfieldPattern}
    pattern.L != transfer.L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    return reciprocaltype(P)( L, pattern.Eθ .* transfer.transfermatrix, pattern.Eϕ .* transfer.transfermatrix )
end

function translate!(incidentfield::P, farfield::F, transfer::PlannedTransfer{C}; reset::Bool=true) where{C, F<:FarfieldPattern, P<:PlaneWaveSpectrum}
    farfield.L != transfer.L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    farfield.L != incidentfield.L && ErrorException("FarFieldPattern and PlaneWaveSpectrum must have the same order L !")
    
    _muladd!(incidentfield.Eθ, farfield.Eθ, transfer.transfermatrix; reset=reset)
    _muladd!(incidentfield.Eϕ, farfield.Eϕ, transfer.transfermatrix; reset=reset)
    return incidentfield
end
# function translate!(incidentfield::P, farfield::F, transfer::OnTheFlyTransfer{L,C}; reset::Bool=true) where{L, C, F<:FarfieldPattern, P<:PlaneWaveSpectrum}  
#     transferplan= PlannedTransfer(transfer.R, transfer.k0, transfer.L)
#     return translate!(incidentfield, farfield, transferplan, reset=reset)
# end
function translate!(incidentfield::PlaneWaveSpectrum, farfield::FarfieldPattern, R::AbstractVector{T}, k₀::T; reset::Bool=true) where{T<:Real}
    transfer=OnTheFlyTransfer{pattern.L, T}(SVector{3,T}(R), k₀)
   return translate!(incidentfield, farfield, transfer, reset=reset)
end


function _adjoint_translate!(incidentfield::P, farfield::F, transfer::PlannedTransfer{C}; reset::Bool=true) where{C, F<:FarfieldPattern, P<:PlaneWaveSpectrum}
    farfield.L != transfer.L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    farfield.L != incidentfield.L && ErrorException("FarFieldPattern and PlaneWaveSpectrum must have the same order L !")
    conj!(transfer.transfermatrix)
    
    _muladd!(farfield.Eθ, incidentfield.Eθ, transfer.transfermatrix; reset=reset)
    _muladd!(farfield.Eϕ, incidentfield.Eϕ, transfer.transfermatrix; reset=reset)

    conj!(transfer.transfermatrix)
    return farfield
end

function _transpose_translate!(incidentfield::P, farfield::F, transfer::PlannedTransfer{C}; reset::Bool=true) where{C, F<:FarfieldPattern, P<:PlaneWaveSpectrum}
    farfield.L != transfer.L && ErrorException("FarFieldPattern and PlannedTransfer must have the same order L !")
    farfield.L != incidentfield.L && ErrorException("FarFieldPattern and PlaneWaveSpectrum must have the same order L !")
    
    _muladd!(farfield.Eθ, incidentfield.Eθ, transfer.transfermatrix; reset=reset)
    _muladd!(farfield.Eϕ, incidentfield.Eϕ, transfer.transfermatrix; reset=reset)
    return incidentfield
end
