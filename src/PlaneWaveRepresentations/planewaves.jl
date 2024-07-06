# const Z₀ = 376.730313669
abstract type PlaneWaveRepresentation <: AntennaFieldRepresentation end

"""
Radiated far-field pattern sampled in angular domain.
Sampling points lie on equiangular grid in ϕ and Gauß-Legendre grid in θ.
Number and location of sampling points is defined by L.
length(θ)=L+1
length(ϕ)=2L+2.

# Fields:
k₀::Float64 : Wavenumber
L::Int64 : Order of the representation
Eθ::Array{<:Complex,2} : matrix storing the complex values for the θ-components
Eϕ::Array{<:Complex,2} : matrix storing the complex values for the ϕ-components
"""
mutable struct FarfieldPattern{C<:Complex} <: PlaneWaveRepresentation
    L::Int64
    Eθ::Array{C,2}
    Eϕ::Array{C,2}
end
import Base.show
function Base.show(io::IO, pattern::FarfieldPattern{C}) where {C}
    L = pattern.L
    print(io, "FarfieldPattern{$C}; L= $L")
end


"""
Incident plane-wave spectrum sampled in angular domain.
Sampling points lie on equiangular grid in ϕ and Gauß-Legendre grid in θ.
Number and location of sampling points is defined by L.
length(θ)=L+1
length(ϕ)=2L+2.

# Fields:
k₀::Float64 : Wavenumber
L::Int64 : Order of the representation
Eθ::Array{Complex{Float64},2} : matrix storing the complex values for the θ-components
Eϕ::Array{Complex{Float64},2} : matrix storing the complex values for the ϕ-components 
"""
mutable struct PlaneWaveSpectrum{C<:Complex} <: PlaneWaveRepresentation
    L::Int64
    Eθ::Array{C,2}
    Eϕ::Array{C,2}
end
import Base.show
function Base.show(io::IO, pattern::PlaneWaveSpectrum{C}) where {C}
    L = pattern.L
    print(io, "PlaneWaveSpectrum{$C}; L= $L")
end

"""
Far-field pattern sampled in angular domain.
Sampling points lie on equiangular grid in ϕ and θ.

# Fields:
k₀::Float64
Jθ::Integer : defines sampling density in θ: Δθ= 2π/Jθ
Jϕ::Integer : defines sampling density in ϕ: Δϕ= 2π/Jϕ
Eθ::Array{Complex{Float64},2} : matrix storing the complex values for the corresponding co polar plane wave
Eϕ::Array{Complex{Float64},2} : matrix storing the complex values for the corresponding cross polar plane wave
"""
mutable struct EquiangularFarFieldPattern{C<:Complex}
    Jθ::Int64
    Jϕ::Int64
    Eθ::Array{C,2} #matrix storing the complex values for the corresponding co polar plane wave
    Eϕ::Array{C,2} #matrix storing the complex values for the corresponding cross polar plane wave
end


"""
Object for a single plane wave

# Fields:
kvec::Array{<:Number,1} : k-vector determining the propagation direction. If complex valued: Evanescent components
pol::Array{<:Number,1} : polarizatiion vector. If complex valued: elliptic polarization
mag::Complex : Complex valued excitation amplitude
"""
mutable struct PlaneWave{C<:Number}
    kvec::Array{Float64,1}
    pol::Array{C,1}
    mag::ComplexF64
end


function converttype(
    T::Type{PlaneWaveSpectrum{C}},
    ff::PlaneWaveRepresentation,
) where {C<:Complex}
    return PlaneWaveSpectrum(ff.L, convert.(C, ff.Eθ), convert.(C, ff.Eϕ))
end
function converttype(
    T::Type{FarfieldPattern{C}},
    pws::PlaneWaveRepresentation,
) where {C<:Complex}
    return FarfieldPattern(pws.L, convert.(C, pws.Eθ), convert.(C, pws.Eϕ))
end

function elementtype(pws::PlaneWaveSpectrum{C}) where {C}
    return C
end
function elementtype(ff::FarfieldPattern{C}) where {C}
    return C
end
function elementtype(pws::PlaneWave{C}) where {C}
    return C
end

function reciprocaltype(T::Type{FarfieldPattern{C}}) where {C}
    return PlaneWaveSpectrum{C}
end
function reciprocaltype(T::Type{PlaneWaveSpectrum{C}}) where {C}
    return FarfieldPattern{C}
end


function ehfield(pw::PlaneWave, R)
    expfac = cis(-udot(R, pw.kvec))
    E = expfac * pw.pol * pw.mag
    H = cross(pw.kvec, E) / Z₀
    return E, H
end
function efield(pw::PlaneWave, R)
    E, _ = ehfield(pw::PlaneWave, R)
    return E
end
function hfield(pw::PlaneWave, R)
    _, H = ehfield(pw::PlaneWave, R)
    return H
end

"""
    revertdirection(pattern::PlaneWaveRepresentation) -> pattern_revert::PlaneWaveRepresentation

Return `PlaneWaveRepresentation` with reverted propagation directions. **F**(**k**) -> **F**(-**k**)
"""
function revertdirection(pattern::P) where {P<:PlaneWaveRepresentation}
    patternout = deepcopy(pattern)
    for k = 1:(pattern.L+1)  # loop over half of ϕ

        # first half of ϕ values
        for kk = 1:(pattern.L+1) # loop over θ
            patternout.Eθ[kk, k] = pattern.Eθ[end-kk+1, k+pattern.L+1]
            patternout.Eϕ[kk, k] = -pattern.Eϕ[end-kk+1, k+pattern.L+1]
        end
        # second half of ϕ values
        for kk = 1:(patternout.L+1) # loop over θ
            patternout.Eθ[kk, k+pattern.L+1] = pattern.Eθ[end-kk+1, k]
            patternout.Eϕ[kk, k+pattern.L+1] = -pattern.Eϕ[end-kk+1, k]
        end
    end
    return patternout
end


"""
    samplingrule(L::Integer) -> ( w::Array{Float64,1}, θ::Array{Float64,1}, ϕ::Array{Float64,1} )

Return sampling points for θ, ϕ, and θ-integration weights w for a far field or plane wave spectrum of order L.
"""
function samplingrule(L::Integer)
    x::Vector{Float64}, w::Vector{Float64} = gausslegendre(L + 1)
    θ = acos.(-x)
    nphi = 2 * L + 2
    ϕ = 2 * pi / (nphi) * collect(0:(nphi-1))

    return w, θ, ϕ
end

function transmission(ff::FarfieldPattern, pws::PlaneWaveSpectrum)
    return transmission(pws, ff)
end
function transmission(ff::FarfieldPattern, pws::PlaneWaveSpectrum, ::Number)
    return transmission(ff, pws)
end
function transmission(pws::PlaneWaveSpectrum, ff::FarfieldPattern)
    if pws.L > ff.L
        return transmission(pws, resample(ff, pws.L))
    elseif pws.L < ff.L
        return transmission(resample(pws, ff.L), ff)
    end

    w, _, _ = samplingrule(ff.L)
    pws_revert = revertdirection(pws)
    return _ewaldintegral(ff, pws_revert, w, ff.L)
    # integrand = (pws_revert.Eθ .* ff.Eθ) + (pws_revert.Eϕ .* ff.Eϕ)
    # w, _, _ = samplingrule(ff.L)

    # return sum(transpose(w) * integrand) * π / (2 * ff.L + 2) / Z₀
end
function transmission(pws::PlaneWaveSpectrum, ff::FarfieldPattern, ::Number)
    return transmission(pws, ff)
end

function _ewaldintegral(
    ff::FarfieldPattern,
    pws::PlaneWaveSpectrum,
    w::Vector{R},
    L::Integer,
) where {R<:Real}
    integrand = (pws.Eθ .* ff.Eθ) + (pws.Eϕ .* ff.Eϕ)
    return sum(transpose(w) * integrand) * π / (2 * L + 2) / Z₀
end

function shiftrepresentation(pattern::PlaneWaveRepresentation, R, k₀)
    newpattern = deepcopy(pattern)
    _, θvec, ϕvec = samplingrule(pattern.L)

    for kϕ in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[kϕ])
        for kθ in eachindex(θvec)
            sint, cost = sincos(θvec[kθ])
            eᵣ = [cosp * sint; sinp * sint; cost]
            ejkr = cis(-k₀ * udot(R, eᵣ))
            newpattern.Eθ[kθ, kϕ] *= ejkr
            newpattern.Eϕ[kθ, kϕ] *= ejkr
        end
    end
    return newpattern
end

"""
    _phaseshiftmatrix!(storagematrix, R, k₀)

In-place calculation of phaseshiftmatrix which needs to be element-wise multiplied to farfield to shift the farfield to new coordinate origin at R.
Assumes that size(storagematrix) mathces the matrices in PlaneWaveRepresentation.
"""
function _phaseshiftmatrix!(storagematrix::AbstractMatrix, R::AbstractVector, k₀::Number)
    a, b = size(storagematrix)
    @assert b == 2 * a
    L = a - 1
    _, θvec, ϕvec = samplingrule(L)
    @assert length(θvec) == a
    @assert length(ϕvec) == b
    Rvec = SVector{3}(R)
    sp = sin.(ϕvec)
    cp = cos.(ϕvec)
    st = sin.(θvec)
    ct = cos.(θvec)
    for kϕ in eachindex(ϕvec)
        sinp, cosp = sp[kϕ], cp[kϕ]
        for kθ in eachindex(θvec)
            sint, cost = st[kθ], ct[kθ]
            eᵣ = typeof(Rvec)(cosp * sint, sinp * sint, cost)
            ejkr = cis(-k₀ * udot(Rvec, eᵣ))
            storagematrix[kθ, kϕ] = ejkr
        end
    end
end


function rotate(
    pattern::PlaneWaveRepresentation,
    χ::Number,
    θ::Number,
    ϕ::Number,
    orderθ::Integer,
    orderϕ::Integer,
)
    R = rot_mat_zyz(-ϕ, -θ, -χ) # rotate sampling points in reverse direction to rotate sampled pattern
    newpattern = deepcopy(pattern)
    _, θvec, ϕvec = samplingrule(pattern.L)

    for kϕ in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[kϕ])
        for kθ in eachindex(θvec)
            sint, cost = sincos(θvec[kθ])

            # calculate unit vectors at sample position
            eᵣ = [cosp * sint; sinp * sint; cost]
            eθ = [cosp * cost; sinp * cost; -sint]
            eϕ = [-sinp; cosp; 0]


            # rotate unit vectors
            eᵣ_rot = R * eᵣ
            eθ_rot = R * eθ
            eϕ_rot = R * eϕ

            # calculate theta and phi angle of rotated sample position
            θ_rot = atan(sqrt(eᵣ_rot[1]^2 + eᵣ_rot[2]^2), eᵣ_rot[3])# result is always positive
            ϕ_rot = mod2pi(atan(eᵣ_rot[2], eᵣ_rot[1]))

            sinpr, cospr = sincos(ϕ_rot)
            sintr, costr = sincos(θ_rot)

            # calculate unit vectors at rotated sample position
            eθ_i = [cospr * costr, sinpr * costr, -sintr]
            eϕ_i = [-sinpr, cospr, 0]

            Eθ, Eϕ = interpolate_single_planewave(
                θ_rot,
                ϕ_rot,
                pattern,
                LocalInterpolation{pattern.L,pattern.L,orderθ,orderϕ,Float64}(),
            )

            newpattern.Eθ[kθ, kϕ] = (eθ_rot ⋅ eϕ_i) * Eϕ + (eθ_rot ⋅ eθ_i) * Eθ
            newpattern.Eϕ[kθ, kϕ] = (eϕ_rot ⋅ eϕ_i) * Eϕ + (eϕ_rot ⋅ eθ_i) * Eθ
        end
    end
    return newpattern
end
function rotate(pattern::PlaneWaveRepresentation, χ::Number, θ::Number, ϕ::Number)
    return rotate(pattern, χ, θ, ϕ, 12, 12)
end
function rotate(
    pattern::PlaneWaveRepresentation;
    χ = 0.0,
    θ = 0.0,
    ϕ = 0.0,
    orderθ = 12,
    orderϕ = 12,
)
    return rotate(pattern, χ, θ, ϕ, orderθ, orderϕ)
end

function efield(p::PlaneWave, R::Array{<:Number})
    return p.mag * p.pol * cis(-udot(R, p.kvec))
end
function hfield(p::PlaneWave, R::Array{<:Number})
    return p.mag / Z₀ * cross(p.kvec / norm(p.kvec), p.pol) * cis(-udot(R, p.kvec))
end
function ehfield(p::PlaneWave, R::Array{<:Number})
    ejkr = cis(-udot(R, p.kvec))
    return p.mag * p.pol * ejkr, p.mag / Z₀ * cross(p.kvec / norm(p.kvec), p.pol) * ejkr
end


function ehfield(
    pws::PlaneWaveSpectrum{C},
    R::A,
    k₀::Number,
) where {A<:AbstractVector{<:Number},C<:Complex}
    T = real(C)
    wθ, θvec, ϕvec = samplingrule(pws.L)
    wϕ = C(0.0, -k₀ / ((4 * pws.L + 4)))

    E, H = SVector{3}(zeros(C, 3)), SVector{3}(zeros(C, 3))
    # H=fill!(Array{Array{C,1}}(undef, size(Rs)), [0,0,0])

    sintcost = sincos.(θvec)
    for kk in eachindex(ϕvec)
        sinp, cosp = convert.(T, sincos(ϕvec[kk]))
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for k in eachindex(θvec)
            sint, cost = convert.(T, sintcost[k])

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Epol = pws.Eθ[k, kk] * eθ + pws.Eϕ[k, kk] * eϕ
            Hpol = pws.Eθ[k, kk] * eϕ - pws.Eϕ[k, kk] * eθ
            ejkr = cis(-k₀ * udot(eᵣ, SVector{3}(R))) * wθ[k] * wϕ
            E += Epol * ejkr
            H += Hpol * ejkr / Z₀
        end
    end
    return E, H
end


function efield(
    pws::PlaneWaveSpectrum{C},
    R::A,
    k₀::Number,
) where {A<:AbstractVector{<:Number},C<:Complex}
    T = real(C)
    wθ, θvec, ϕvec = samplingrule(pws.L)
    wϕ = C(0.0, -k₀ / ((4 * pws.L + 4)))

    E = SVector{3}(zeros(C, 3))
    # H=fill!(Array{Array{C,1}}(undef, size(Rs)), [0,0,0])

    sintcost = sincos.(θvec)
    for kk in eachindex(ϕvec)
        sinp, cosp = convert.(T, sincos(ϕvec[kk]))
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for k in eachindex(θvec)
            sint, cost = convert.(T, sintcost[k])

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Epol = pws.Eθ[k, kk] * eθ + pws.Eϕ[k, kk] * eϕ
            ejkr = cis(-k₀ * udot(eᵣ, SVector{3}(R))) * wθ[k] * wϕ
            E += Epol * ejkr
        end
    end
    return E
end


function hfield(
    pws::PlaneWaveSpectrum{C},
    R::A,
    k₀::Number,
) where {A<:AbstractVector{<:Number},C<:Complex}
    T = real(C)
    wθ, θvec, ϕvec = samplingrule(pws.L)
    wϕ = C(0.0, -k₀ / ((4 * pws.L + 4)))

    H = SVector{3}(zeros(C, 3))
    # H=fill!(Array{Array{C,1}}(undef, size(Rs)), [0,0,0])
    sintcost = sincos.(θvec)
    for kk in eachindex(ϕvec)
        sinp, cosp = convert.(T, sincos(ϕvec[kk]))
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for k in eachindex(θvec)
            sint, cost = convert.(T, sintcost[k])

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Hpol = pws.Eθ[k, kk] * eϕ - pws.Eϕ[k, kk] * eθ
            ejkr = cis(-k₀ * udot(eᵣ, SVector{3}(R))) * wθ[k] * wϕ
            H += Hpol * ejkr / Z₀
        end
    end
    return H
end


# function efield(pws::PlaneWaveSpectrum, R::A, k₀::Number) where {A<:AbstractVector{<:Number}}
#     fieldtuples = ehfield(pws, R, k₀)
#     return fieldtuples[1]
# end


# function hfield(pws::PlaneWaveSpectrum, R::A, k₀::Number) where {A<:AbstractVector{<:Number}}
#     fieldtuples = ehfield(pws, R, k₀)
#     return fieldtuples[2]
# end



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
        floatℓ = T(ℓ)
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

#TODO: Better name: _add_or_overwrite!
function _add!(
    storage::PlaneWaveRepresentation,
    summand::PlaneWaveRepresentation;
    reset::Bool = false,
)
    if reset
        storage.Eθ .= summand.Eθ
        storage.Eϕ .= summand.Eϕ
    else
        storage.Eθ .+= summand.Eθ
        storage.Eϕ .+= summand.Eϕ
    end
    return storage
end
function _add!(storage::AbstractMatrix, summand::AbstractMatrix; reset::Bool = false)
    if reset
        storage .= summand
    else
        storage .+= summand
    end
    return storage
end

#TODO: Better name: _muladd_or_muloverwrite!
function _muladd!(
    storage::AbstractMatrix,
    summand::AbstractMatrix,
    factor::AbstractMatrix;
    reset::Bool = false,
)
    if reset
        storage .= factor .* summand
    else
        storage .= storage .+ (factor .* summand)
    end
    return storage
end
function _muladd!(
    storage::PlaneWaveRepresentation,
    summand::PlaneWaveRepresentation,
    factor::Number;
    reset::Bool = false,
)
    if reset
        storage.Eθ .= factor .* summand.Eθ
        storage.Eϕ .= factor .* summand.Eϕ
    else
        storage.Eθ .+= factor .* summand.Eθ
        storage.Eϕ .+= factor .* summand.Eϕ
    end
    return storage
end

function _mul!(
    storage::PlaneWaveRepresentation,
    factor::PlaneWaveRepresentation;
    reset::Bool = false,
)
    if reset
        storage.Eθ .= factor.Eθ
        storage.Eϕ .= factor.Eϕ
    else
        storage.Eθ .*= factor.Eθ
        storage.Eϕ .*= factor.Eϕ
    end
    return storage
end

function _mul!(storage::PlaneWaveRepresentation, factor::Number; reset::Bool = false)
    if reset
        fill!(storage.Eθ, factor)
        fill!(storage.Eθ, factor)
    else
        storage.Eθ .*= factor
        storage.Eϕ .*= factor
    end
    return storage
end
function _mul!(factor::Number, storage::PlaneWaveRepresentation; reset::Bool = false)
    _mul!(storage, factor, reset = reset)
    return storage
end
