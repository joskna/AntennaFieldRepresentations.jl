
"""
    PlaneWaveExpansion{P,Y,C} <: AntennaFieldRepresentation{P,C}

Collection of electromagnetic plane waves propagating into various directions.

Behaves like an `AbstractVector{C}` with extra context.

# Type Parameters
- `P <: PropagationType`
- `Y <: SphereSamplingStrategy` : Defines the propagation directions of the  stored plane wave samples
- `C <: Number`
"""
struct PlaneWaveExpansion{P<:PropagationType,Y<:SphereSamplingStrategy,C} <:
       AntennaFieldRepresentation{P,C}
    samplingstrategy::Y
    EθEϕ::Array{C,3}
    wavenumber::Number
    buffer::Vector{C}
end

"""
    PlaneWaveExpansion(P, samplingstrategy, Eθ, Eϕ, wavenumber)

Returns a collection of electromagnetic plane waves propagating into various directions.

# Arguments:
- `P::PropagationType`: Radiated(), Incident(), or Absorbed()
- `samplingstrategy<:SphereSamplingStrategy`: sampling strategy which defines the propagation directions of the plane wave samples
- `Eθ::Matrix{C}`: θ-component amplitudes of the plane waves
- `Eϕ::Matrix{C}`: ϕ-component amplitudes of the plane waves
- `wavenumber`: wavenumber
"""
function PlaneWaveExpansion(
    propagation::P,
    samplingstrategy::S,
    Eθ::Matrix{C},
    Eϕ::Matrix{C},
    wavenumber::Number,
) where {P<:PropagationType,S<:SphereSamplingStrategy,C}
    a, b = size(Eθ)
    a_, b_ = size(Eϕ)
    !(a == a_ && b == b_) && error("Dimensions of Eθ and Eϕ must be equal")
    Eθϕ = Array{C}(undef, a, b, 2)
    Eθϕ[:, :, 1] .= Eθ
    Eθϕ[:, :, 2] .= Eϕ
    buffer = reshape(Eθϕ, length(Eθϕ))
    return PlaneWaveExpansion{P,S,C}(samplingstrategy, Eθϕ, wavenumber, buffer)
end
function PlaneWaveExpansion(
    ::P,
    samplingstrategy::S,
    wavenumber::Number,
) where {P<:PropagationType,S<:SphereSamplingStrategy}
    a, b = _countsamples(samplingstrategy)
    Eθϕ = Array{ComplexF64}(undef, a, b, 2)
    buffer = reshape(Eθϕ, length(Eθϕ))
    return PlaneWaveExpansion{P,S,ComplexF64}(samplingstrategy, Eθϕ, wavenumber, buffer)
end
# function _eθ(p::PlaneWaveExpansion)
#     _, s2 = _count_samples(p.samplingstrategy)
#     return view(p.EθEϕ, :, 1:s2)
# end
function _eθ(p::PlaneWaveExpansion{P,Y,C}) where {P,Y,C}

    return view(p.EθEϕ, :, :, 1)
end
# function _eϕ(p::PlaneWaveExpansion)
#     _, s2 = _count_samples(p.samplingstrategy)
#     return view(p.EθEϕ, :, (s2+1):2*s2)
# end
function _eϕ(p::PlaneWaveExpansion)
    return view(p.EθEϕ, :, :, 2)
end
function asvector(p::PlaneWaveExpansion)
    return p.buffer
end
Base.size(p::PlaneWaveExpansion) = size(asvector(p))
Base.getindex(p::PlaneWaveExpansion, i) = Base.getindex(asvector(p), i)
Base.setindex!(p::PlaneWaveExpansion, i, v) = Base.setindex!(asvector(p), i, v)
function Base.similar(p::PlaneWaveExpansion{P,S,C}) where {P,S,C}
    EθEϕ = similar(p.EθEϕ)
    return PlaneWaveExpansion{P,S,C}(
        p.samplingstrategy,
        EθEϕ,
        p.wavenumber,
        reshape(EθEϕ, length(EθEϕ)),
    )
end
function setwavenumber!(p::PlaneWaveExpansion{P,S,C}, val) where {P,C,S}
    swe = PlaneWaveExpansion{P,C,S}(p.samplingstrategy, p.EθEϕ, val, p.buffer)
    return swe
end

function efield!(
    storage,
    pwe::PlaneWaveExpansion{Incident,S,C},
    R;
    reset = true,
) where {S<:SphereSamplingStrategy,C<:Complex}
    θweights, ϕweights, θs, ϕs = weightsandsamples(pwe.samplingstrategy)

    T = real(C)
    Rvec = SVector{3}(R)
    reset && fill!(storage, zero(C))
    sintcost = sincos.(θs)
    for (kk, ϕ) in enumerate(ϕs)
        sinp, cosp = sincos(ϕ)
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for (k, θ) in enumerate(θs)
            sint, cost = sintcost[k]

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Epol = _eθ(pwe)[k, kk] * eθ + _eϕ(pwe)[k, kk] * eϕ
            ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
            storage += Epol * ejkr
        end
    end
    return storage
end

function hfield!(
    storage,
    pwe::PlaneWaveExpansion{Incident,S,C},
    R;
    reset = true,
) where {S<:SphereSamplingStrategy,C<:Complex}
    θweights, ϕweights, θs, ϕs = weightsandsamples(pwe.samplingstrategy)

    T = real(C)
    Rvec = SVector{3}(R)
    reset && fill!(storage, zero(C))
    sintcost = sincos.(θs)
    for (kk, ϕ) in enumerate(ϕs)
        sinp, cosp = sincos(ϕ)
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for (k, θ) in enumerate(θs)
            sint, cost = sintcost[k]

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Hpol = _eθ(pwe)[k, kk] * eϕ - _eϕ(pwe)[k, kk] * eθ
            ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
            storage += Hpol * ejkr / Z₀
        end
    end
    return storage
end

function ehfield!(
    estorage,
    hstorage,
    pwe::PlaneWaveExpansion{Incident,S,C},
    R;
    reset = true,
) where {S<:SphereSamplingStrategy,C<:Complex}
    θweights, ϕweights, θs, ϕs = weightsandsamples(pwe.samplingstrategy)

    T = real(C)
    Rvec = SVector{3}(R)
    reset && fill!(estorage, zero(C))
    reset && fill!(hstorage, zero(C))
    sintcost = sincos.(θs)
    for (kk, ϕ) in enumerate(ϕs)
        sinp, cosp = sincos(ϕ)
        eϕ = SVector{3}(-sinp, cosp, zero(T))
        for (k, θ) in enumerate(θs)
            sint, cost = sintcost[k]

            eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
            eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

            Hpol = _eθ(pwe)[k, kk] * eϕ - _eϕ(pwe)[k, kk] * eθ
            Epol = _eθ(pwe)[k, kk] * eθ + _eϕ(pwe)[k, kk] * eϕ
            ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
            storage .+= Epol * ejkr
            hstorage .+= Hpol * ejkr / Z₀
        end
    end
    return estorage, hstorage
end

function _add_or_reset!(
    storage::PlaneWaveExpansion,
    summand::PlaneWaveExpansion;
    reset::Bool = false,
)
    if reset
        storage.EθEϕ .= summand.EθEϕ
    else
        storage.EθEϕ .+= summand.EθEϕ
    end
    return storage
end
function _add_or_reset!(
    storage::AbstractMatrix,
    summand::AbstractMatrix;
    reset::Bool = false,
)
    if reset
        storage .= summand
    else
        storage .+= summand
    end
    return storage
end


function _muladd_or_mulreset!(
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
function _muladd_or_mulreset!(
    storage::PlaneWaveExpansion,
    summand::PlaneWaveExpansion,
    factor::Number;
    reset::Bool = false,
)
    if reset
        storage.EθEϕ .= summand.EθEϕ
        storage.EθEϕ .*= factor
    else
        storage.EθEϕ .+= factor .* summand.EθEϕ
    end
    return storage
end

function _muladd_or_mulreset!(
    storage::AbstractVector,
    operator,
    summand::AbstractVector;
    reset::Bool = false,
)
    if reset
        storage .= operator * summand
    else
        storage .+= operator * summand
    end
    return storage
end

function _mul_or_reset!(
    storage::PlaneWaveExpansion,
    factor::PlaneWaveExpansion;
    reset::Bool = false,
)
    if reset
        storage.EθEϕ .= factor.EθEϕ
    else
        storage.EθEϕ .*= factor.EθEϕ
    end
    return storage
end

function _mul_or_reset!(storage::PlaneWaveExpansion, factor::Number; reset::Bool = false)
    if reset
        fill!(storage.EθEϕ, factor)
    else
        storage.EθEϕ .*= factor
    end
    return storage
end
function _mul_or_reset!(factor::Number, storage::PlaneWaveExpansion; reset::Bool = false)
    _mul_or_reset!(storage, factor, reset = reset)
    return storage
end

function equivalentorder(pwe::PlaneWaveExpansion; ϵ = 1e-7)
    a, b, c = size(pwe.EθEϕ)
    L = minimum([a - 1, (b - 1) ÷ 2])
    return L
end

"""
    _phaseshiftmatrix!(storagematrix, R, k₀)

In-place calculation of phaseshiftmatrix which needs to be element-wise multiplied to farfield to shift the farfield to new coordinate origin at R.
Assumes that size(storagematrix) mathces the matrices in PlaneWaveRepresentation.
"""
function _phaseshiftmatrix!(
    storagematrix::AbstractMatrix,
    R,
    k₀::Number,
    sampling::Y,
) where {Y<:SphereSamplingStrategy}
    T = Float64
    a, b = size(storagematrix)
    θvec, ϕvec = samples(sampling)
    L = a - 1
    @assert length(θvec) == a
    @assert length(ϕvec) == b

    sp = T.(sin.(ϕvec))
    cp = T.(cos.(ϕvec))
    st = T.(sin.(θvec))
    ct = T.(cos.(θvec))
    eᵣ = zeros(T, 3)
    for kϕ in eachindex(ϕvec)
        sinp, cosp = sp[kϕ], cp[kϕ]
        for kθ in eachindex(θvec)
            sint, cost = st[kθ], ct[kθ]
            eᵣ .= cosp * sint, sinp * sint, cost
            ejkr = cis(-k₀ * udot(R, eᵣ))
            storagematrix[kθ, kϕ] = ejkr
        end
    end
    return storagematrix
end

include("interpolation.jl")

function rotate(
    pattern::W,
    χ::T,
    θ::T,
    ϕ::T;
    orderθ::Integer = 12,
    orderϕ::Integer = 12,
) where {W<:PlaneWaveExpansion,T<:Number}
    return rotate!(similar(pattern), pattern, χ, θ, ϕ; orderθ = orderθ, orderϕ = orderϕ)
end
function rotate!(
    rotated_pattern::W,
    pattern::W,
    χ::T,
    θ::T,
    ϕ::T;
    orderθ::Integer = 12,
    orderϕ::Integer = 12,
) where {W<:PlaneWaveExpansion,T<:Number}
    R = rot_mat_zyz(-ϕ, -θ, -χ) # rotate sampling points in reverse direction to rotate sampled pattern
    θvec, ϕvec = samples(pattern.samplingstrategy)

    for kϕ in eachindex(ϕvec)
        sinp, cosp = sincos(ϕvec[kϕ])
        for kθ in eachindex(θvec)
            sint, cost = sincos(θvec[kθ])

            # calculate unit vectors at sample position
            eᵣ = SVector{3}([cosp * sint; sinp * sint; cost])
            eθ = SVector{3}([cosp * cost; sinp * cost; -sint])
            eϕ = SVector{3}([-sinp; cosp; 0])


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
                (θ_rot, ϕ_rot),
                pattern;
                orderθ = orderθ,
                orderϕ = orderϕ,
            )

            rotated_pattern.EθEϕ[kθ, kϕ, 1] = (eθ_rot ⋅ eϕ_i) * Eϕ + (eθ_rot ⋅ eθ_i) * Eθ
            rotated_pattern.EθEϕ[kθ, kϕ, 2] = (eϕ_rot ⋅ eϕ_i) * Eϕ + (eϕ_rot ⋅ eθ_i) * Eθ
        end
    end
    return rotated_pattern
end

function interpolate(
    θϕ::Tuple{T,T},
    pattern::PlaneWaveExpansion;
    orderθ = 12,
    orderϕ = 12,
) where {T}
    Eθ, Eϕ = interpolate_single_planewave(θϕ, pattern, orderθ = orderθ, orderϕ = orderϕ)
    return Eθ, Eϕ
end

include("resampling.jl")
function resample(
    newsamplingstrategy::Y1,
    pwe::PlaneWaveExpansion{P,Y2,C},
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,P<:PropagationType,C}
    rsm = ResampleMap(newsamplingstrategy, pwe.samplingstrategy)
    y = rsm * pwe

    θs, ϕs = samples(newsamplingstrategy)
    EθEϕ = reshape(y, length(θs), length(ϕs), 2)
    return PlaneWaveExpansion{P,Y1,C}(newsamplingstrategy, EθEϕ, pwe.wavenumber, y)
end

include("transfer.jl")
