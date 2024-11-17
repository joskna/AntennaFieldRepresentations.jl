"""
    AbstractSphericalCoefficients{C <: Number }

Supertype for collections of spherical coefficients.

Can be accessed like an AbstractVector{C} with a single index `j` or with a triple index `(s,ℓ,m)`, where `j = 2 * (ℓ * (ℓ + 1) + m - 1) + s`.
"""
abstract type AbstractSphericalCoefficients{C<:Number} <: AbstractVector{C} end

"""
    SphericalCoefficients{C} <: AbstractSphericalCoefficients{C}

Collections of coefficients for a spherical wave expansion. 
    
Can be accessed like an AbstractVector{C} with a single index `j` or with a triple index `(s,ℓ,m)`, where `j = 2 * (ℓ * (ℓ + 1) + m - 1) + s`.

# Examples
```julia
values=ComplexF64.(collect(1:16))
sph_coefficients= SphericalCoefficients(values)
println(sph_coefficients[9])
println(sph_coefficients[1, 2, -1])
sph_coefficients[1, 2, -1] = 111.0im
println(sph_coefficients[9])
```
"""
struct SphericalCoefficients{C} <: AbstractSphericalCoefficients{C}
    coefficients::Vector{C}
end
function asvector(swe::SphericalCoefficients)
    return swe.coefficients
end
Base.getindex(swe::SphericalCoefficients, j) = getindex(swe.coefficients, j)
Base.getindex(swe::SphericalCoefficients, s, ℓ, m) =
    getindex(swe.coefficients, sℓm_to_j(s, ℓ, m))
Base.setindex!(swe::SphericalCoefficients, v, j) = setindex!(swe.coefficients, v, j)
Base.setindex!(swe::SphericalCoefficients, v, s, ℓ, m) =
    setindex!(swe.coefficients, v, sℓm_to_j(s, ℓ, m))
function Base.similar(swe::SphericalCoefficients)
    return SphericalCoefficients(similar(swe.coefficients))
end
function Base.size(swe::SphericalCoefficients)
    return size(swe.coefficients)
end
function SphericalCoefficients(coefficients::AbstractVector{C}) where {C<:Number}
    return SphericalCoefficients{C}(convert(Vector{C}, coefficients))
end

"""
    sℓm_to_j(s,ℓ,m)

    Convert multi-index s ℓ m to single index j in spherical wave expansion
"""
function sℓm_to_j(s, ℓ, m)
    j = 2 * (ℓ * (ℓ + 1) + m - 1) + s
    return j

end

"""
    FirstOrderSphericalCoefficients{C} <: AbstractSphericalCoefficients{C}

Sparse representation of a collections of only first-order coefficients for a spherical wave expansion. 
    
Can be accessed like an AbstractVector{C} with a single index `j` or with a triple index `(s,ℓ,m)`, where `j = 2 * (ℓ * (ℓ + 1) + m - 1) + s`.
All coefficients with `m ≠ ± 1` are zero by definition.

# Examples
```julia
values=ComplexF64.(collect(1:16))
sph_coefficients= SphericalCoefficients(values)
firstorder_sph_coefficients = FirstOrderSphericalCoefficients(sph_coefficients)
println(firstorder_sph_coefficients[1, 2, -1])
firstorder_sph_coefficients[1, 2, -1] = 111.0 im 
println(firstorder_sph_coefficients[9])
```
"""
struct FirstOrderSphericalCoefficients{C} <: AbstractSphericalCoefficients{C}
    coefficients1ℓplus::Vector{C}
    coefficients1ℓminus::Vector{C}
    coefficients2ℓplus::Vector{C}
    coefficients2ℓminus::Vector{C}
end
function Base.size(swe::FirstOrderSphericalCoefficients)
    lmax = length(swe.coefficients1ℓminus)
    return (sℓm_to_j(2, lmax, lmax),)
end
function Base.getindex(swe::FirstOrderSphericalCoefficients{C}, j) where {C}
    s, ℓ, m = j_to_sℓm(j)
    return getindex(swe, s, ℓ, m)
end
function Base.getindex(swe::FirstOrderSphericalCoefficients{C}, s, ℓ, m) where {C}
    abs(m) != 1 && return zero(C)
    s != 1 && s != 2 && error("Index s must be 1 or 2")
    abs(m) > ℓ && error("|m| must be smaller or equal to ℓ.")
    s == 1 && m == 1 && return swe.coefficients1ℓplus[ℓ]
    s == 2 && m == 1 && return swe.coefficients2ℓplus[ℓ]
    s == 1 && m == -1 && return swe.coefficients1ℓminus[ℓ]
    s == 2 && m == -1 && return swe.coefficients2ℓminus[ℓ]
end
function Base.setindex!(swe::FirstOrderSphericalCoefficients, v, j)
    s, ℓ, m = j_to_sℓm(j)
    return setindex!(swe, v, s, ℓ, m)
end
function Base.setindex!(swe::FirstOrderSphericalCoefficients, v, s, ℓ, m)
    #    abs(m) != 1 && @warn "Setting index with m ≠ ±1  is invalid in a first order expansion. Index was not set. Perhaps you want to convert to non-first-order `SphericalCoefficients`?"
    s != 1 && s != 2 && error("Index s must be 1 or 2")
    abs(m) > ℓ && error("|m| must be smaller or equal to ℓ.")
    s == 1 && m == 1 && setindex!(swe.coefficients1ℓplus, v, ℓ)
    s == 2 && m == 1 && setindex!(swe.coefficients2ℓplus, v, ℓ)
    s == 1 && m == -1 && setindex!(swe.coefficients1ℓminus, v, ℓ)
    s == 2 && m == -1 && setindex!(swe.coefficients2ℓminus, v, ℓ)
end
function Base.similar(swe::FirstOrderSphericalCoefficients{C}) where {C}
    return FirstOrderSphericalCoefficients{C}(
        similar(swe.coefficients1ℓplus),
        similar(swe.coefficients1ℓminus),
        similar(swe.coefficients2ℓplus),
        similar(swe.coefficients2ℓminus),
    )
end
function asvector(swe::FirstOrderSphericalCoefficients{C}) where {C}
    len = length(swe)
    coefficients = zeros(C, len)
    _, Lmax, __ = j_to_sℓm(len)
    for ℓ = 1:Lmax
        firstindex = sℓm_to_j(1, ℓ, 1)
        coefficients[firstindex] = swe.coefficients1ℓminus[ℓ]
        coefficients[firstindex+1] = swe.coefficients2ℓminus[ℓ]
        coefficients[firstindex+5] = swe.coefficients1ℓplus[ℓ]
        coefficients[firstindex+6] = swe.coefficients2ℓplus[ℓ]
    end
    return coefficients
end

"""
    FirstOrderSphericalCoefficients(coefficients::SphericalCoefficients)

Return sparse representation of a collections of only first-order spherical coefficients.
All non-first-order coefficients (index m ≠ ±1) are dropped and treated as zero.

# Examples
```julia
values=ComplexF64.(collect(1:16))
sph_coefficients= SphericalCoefficients(values)
firstorder_sph_coefficients = FirstOrderSphericalCoefficients(sph_coefficients)
println(firstorder_sph_coefficients[1, 2, -1])
println(firstorder_sph_coefficients[1, 2, -2])
```
"""
function FirstOrderSphericalCoefficients(swe::SphericalCoefficients{C}) where {C<:Number}
    _, L, __ = j_to_sℓm(length(swe.coefficients))
    coefficients1ℓplus = zeros(C, L)
    coefficients1ℓminus = zeros(C, L)
    coefficients2ℓplus = zeros(C, L)
    coefficients2ℓminus = zeros(C, L)
    for ℓ = 1:L
        coefficients1ℓplus[ℓ] = swe[1, ℓ, 1]
        coefficients1ℓminus[ℓ] = swe[1, ℓ, -1]
        coefficients2ℓplus[ℓ] = swe[2, ℓ, 1]
        coefficients2ℓminus[ℓ] = swe[2, ℓ, -1]
    end
    return FirstOrderSphericalCoefficients{C}(
        coefficients1ℓplus,
        coefficients1ℓminus,
        coefficients2ℓplus,
        coefficients2ℓminus,
    )
end
function FirstOrderSphericalCoefficients(coefficients::AbstractVector{C}) where {C<:Number}
    _, L, __ = j_to_sℓm(length(coefficients))
    coefficients1ℓplus = zeros(C, L)
    coefficients1ℓminus = zeros(C, L)
    coefficients2ℓplus = zeros(C, L)
    coefficients2ℓminus = zeros(C, L)
    for ℓ = 1:L
        coefficients1ℓplus[ℓ] = coefficients[sℓm_to_j(1, ℓ, 1)]
        coefficients1ℓminus[ℓ] = coefficients[sℓm_to_j(1, ℓ, -1)]
        coefficients2ℓplus[ℓ] = coefficients[sℓm_to_j(2, ℓ, 1)]
        coefficients2ℓminus[ℓ] = coefficients[sℓm_to_j(2, ℓ, -1)]
    end
    return FirstOrderSphericalCoefficients{C}(
        coefficients1ℓplus,
        coefficients1ℓminus,
        coefficients2ℓplus,
        coefficients2ℓminus,
    )
end

"""
    j_to_sℓm(j)

    Convert single index j to multi-index s ℓ m in spherical wave expansion
"""
function j_to_sℓm(j::Integer)
    jtype = typeof(j)
    s = 0
    if isodd(j)
        s = 1
    else
        s = 2
    end
    ℓ = floor(sqrt((j - s) / 2 + 1))
    m = (j - s) / 2 + 1 - ℓ * (ℓ + 1)
    return jtype(s), jtype(ℓ), jtype(m)
end

"""
    SphericalWaveExpansion{P,H,C} <: AntennaFieldRepresentation{P, C}

Representation of an electromagnetic field as superposition spherical vector wave functions.

Behaves like an `AbstractVector{C}` with extra context.

# Type parameters
- `P <: PropagationType`
- `H <: AbstractSphericalCoefficients{C <: Number}` : defines how the spherical coefficients are ordered in memory.
- `C <: Complex`
"""
struct SphericalWaveExpansion{
    P<:PropagationType,
    H<:AbstractSphericalCoefficients,
    C<:Number,
} <: AntennaFieldRepresentation{P,C}
    coefficients::H
    wavenumber::Number
end
function asvector(s::SphericalWaveExpansion)
    return asvector(s.coefficients)
end
Base.getindex(swe::SphericalWaveExpansion, j) = getindex(swe.coefficients, j)
Base.getindex(swe::SphericalWaveExpansion, s, ℓ, m) =
    getindex(swe.coefficients, sℓm_to_j(s, ℓ, m))
Base.setindex!(swe::SphericalWaveExpansion, v, j) = setindex!(swe.coefficients, v, j)
Base.setindex!(swe::SphericalWaveExpansion, v, s, ℓ, m) =
    setindex!(swe.coefficients, v, sℓm_to_j(s, ℓ, m))
Base.size(swe::SphericalWaveExpansion) = size(swe.coefficients)
function Base.similar(swe::SphericalWaveExpansion{P,H,C}) where {P,H,C}
    return SphericalWaveExpansion{P,H,C}(similar(swe.coefficients), swe.wavenumber)
end
function SphericalWaveExpansion(
    ::P,
    coefficients::AbstractVector{C},
    wavenumber::Number,
) where {P<:PropagationType,C}
    H = SphericalCoefficients{C}
    return SphericalWaveExpansion{P,H,C}(H(coefficients), wavenumber)
end
function SphericalWaveExpansion(
    ::P,
    coefficients::H,
    wavenumber::Number,
) where {P<:PropagationType,C,H<:AbstractSphericalCoefficients{C}}
    return SphericalWaveExpansion{P,H,C}(coefficients, wavenumber)
end
function setwavenumber!(swe::SphericalWaveExpansion{P,H,C}, val) where {P,H,C}
    swe = SphericalWaveExpansion{P,H,C}(swe.coefficients, val)
    return swe
end

include("normalizedlegendre.jl")
include("radialfunctions.jl")
include("modefunctions.jl")

include("fastspherical.jl")

function efield!(
    storage,
    aut_field::SphericalWaveExpansion{P,H,C},
    R;
    reset = true,
) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients{C}}
    sqrtZ₀ = convert(C, sqrt(Z₀))
    k0 = getwavenumber(aut_field)
    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(aut_field)

        Fx, Fy, Fz = F_sℓm_cartesian_array(J, P(), R, k0)
        fac = C(k0) * sqrtZ₀
        # Ex= fac* udot(Fx, aut_field)
        # Ey= fac* udot(Fy, aut_field) 
        # Ez= fac* udot(Fz, aut_field)  
        # E= [Ex, Ey, Ez]
        E = fac .* (udot(F, aut_field) for F in (Fx, Fy, Fz))

    else
        E = _E_at_origin(aut_field)
    end

    if reset
        storage .= E
    else
        storage .= E .+ storage
    end

    return storage

end

function hfield!(
    storage,
    aut_field::SphericalWaveExpansion{P,H,C},
    R;
    reset = true,
) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients{C}}
    sqrtZ₀ = convert(C, sqrt(Z₀))
    k0 = getwavenumber(aut_field)
    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(aut_field)

        Fx, Fy, Fz = curlF_sℓm_cartesian_array(J, P(), R, k0)
        fac = (C(0.0, k0) / sqrtZ₀)
        Hfield = fac .* (udot(F, aut_field) for F in (Fx, Fy, Fz))

    else
        Hfield = _H_at_origin(aut_field)
    end

    if reset
        storage .= Hfield
    else
        storage .= Hfield .+ storage
    end
    return storage
end

function ehfield!(
    storage_efield,
    storage_hfield,
    aut_field::SphericalWaveExpansion{P,H,C},
    R;
    reset = true,
) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients{C}}
    sqrtZ₀ = C(sqrt(Z₀))
    k0 = getwavenumber(aut_field)
    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(aut_field)

        Fx, Fy, Fz = F_sℓm_cartesian_array(J, P(), R, k0)
        Ecartesian = (C(k0) * sqrtZ₀) .* (udot(F, aut_field) for F in (Fx, Fy, Fz))

        Hx = udot(Fx[1:2:J], aut_field[2:2:J]) + udot(Fx[2:2:J], aut_field[1:2:J])
        Hy = udot(Fy[1:2:J], aut_field[2:2:J]) + udot(Fy[2:2:J], aut_field[1:2:J])
        Hz = udot(Fz[1:2:J], aut_field[2:2:J]) + udot(Fz[2:2:J], aut_field[1:2:J])
        Hcartesian = (C(0.0, k0) / sqrtZ₀) .* [Hx, Hy, Hz]

    else
        Ecartesian = E_at_origin(aut_field)
        Hcartesian = H_at_origin(aut_field)
    end

    if reset
        storage_efield .= Ecartesian
        storage_hfield .= Hcartesian
    else
        storage_efield .= Ecartesian .+ storage_efield
        storage_hfield .= Hcartesian .+ storage_hfield
    end
    return storage_efield, storage_hfield
end

function _H_at_origin(
    aut_field::SphericalWaveExpansion{Incident,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    F2m11, F201, F211 = _F2m1cartesian_at_origin(C)
    Hcartesian =
        C(0, getwavenumber(aut_field)) / sqrt(Z₀) * (
            F211 * aut_field[1, 1, 1] +
            F2m11 * aut_field[1, 1, -1] +
            F201 * aut_field[1, 1, 0]
        )
    return Hcartesian
end
function _H_at_origin(
    aut_field::SphericalWaveExpansion{Absorbed,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    return _H_at_origin(
        SphericalWaveExpansion{Incident,H,C}(aut_field.coefficients, aut_field.wavenumber),
    ) .+ C(0.0, -Inf)
end
function _H_at_origin(
    aut_field::SphericalWaveExpansion{Radiated,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    return _H_at_origin(
        SphericalWaveExpansion{Incident,H,C}(aut_field.coefficients, aut_field.wavenumber),
    ) .+ C(0.0, Inf)
end

function _E_at_origin(
    aut_field::SphericalWaveExpansion{Incident,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    F2m11, F201, F211 = _F2m1cartesian_at_origin(C)
    Ecartesian =
        getwavenumber(aut_field) *
        sqrt(Z₀) *
        (
            F211 * aut_field[2, 1, 1] +
            F2m11 * aut_field[2, 1, -1] +
            F201 * aut_field[2, 1, 0]
        )
    return Ecartesian
end
function _E_at_origin(
    aut_field::SphericalWaveExpansion{Absorbed,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    return _E_at_origin(
        SphericalWaveExpansion{Incident,H,C}(aut_field.coefficients, aut_field.wavenumber),
    ) .+ C(0.0, -Inf)
end
function _E_at_origin(
    aut_field::SphericalWaveExpansion{Radiated,H,C},
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    return _E_at_origin(
        SphericalWaveExpansion{Incident,H,C}(aut_field.coefficients, aut_field.wavenumber),
    ) .+ C(0.0, Inf)
end


function farfield(
    aut_field::SphericalWaveExpansion{Radiated,H,C},
    ϑ::Number,
    φ::Real,
) where {C<:Number,H<:AbstractSphericalCoefficients{C}}
    J = length(aut_field)
    Kϑ, Kφ = K_sℓm_array(J, ϑ, φ)
    Fϑ = udot(Kϑ[1:J], aut_field)
    Fφ = udot(Kφ[1:J], aut_field)
    return sqrt(Z₀) * Fϑ, sqrt(Z₀) * Fφ
end


function rotate!(
    rotated_aut_field::SphericalWaveExpansion{P,H,C},
    aut_field::SphericalWaveExpansion{P,H,C},
    χ::Real,
    θ::Real,
    ϕ::Real,
) where {P<:PropagationType,C<:Number,H<:AbstractSphericalCoefficients{C}}
    _, Lmax, __ = j_to_sℓm(length(aut_field))

    fill!(rotated_aut_field.coefficients, 0.0)

    for ℓ = 1:Lmax
        d = wignerd(ℓ, -θ)
        for m = (-ℓ):ℓ
            expϕ = cis(-m * ϕ)
            for s = 1:2
                jj = sℓm_to_j(s, ℓ, m)
                for μ = (-ℓ):ℓ
                    j = sℓm_to_j(s, ℓ, μ)
                    if j <= length(aut_field)

                        rotated_aut_field.coefficients[jj] +=
                            expϕ * d[μ+ℓ+1, m+ℓ+1] * cis(-μ * χ) * aut_field.coefficients[j]
                    end
                end
            end
        end
    end

    return rotated_aut_field
end




function αtoβ(α::AbstractVector)
    β = similar(α, length(α))
    return αtoβ!(β, α)
end
function αtoβ!(β, α::AbstractVector)
    for k = 1:length(α)
        s, ℓ, m = j_to_sℓm(k)
        # β[k] = (-1)^(m) * (α[sℓm_to_j(s, ℓ, -m)])
        β[k] = _negpow1(m) * (α[sℓm_to_j(s, ℓ, -m)])
    end
    β .*= 0.5
    return β
end

function βtoα!(α, β::AbstractVector)
    for k = 1:length(β)
        s, ℓ, m = j_to_sℓm(k)
        # α[k] = (-1)^(m) * (β[sℓm_to_j(s, ℓ, -m)])
        α[k] = _negpow1(m) * (β[sℓm_to_j(s, ℓ, -m)])
    end
    α .= α .* 2
    return α
end
function βtoα(β::AbstractVector)
    α = similar(β, length(β))
    return βtoα!(α, β)
end

function _negpow1(m::Integer)
    if isodd(m)
        return -1
    else
        return 1
    end
end



function _outputmode_dipo2sph(Psph::PropagationType, Pdip::PropagationType)
    if Pdip == Radiated()
        return _dualtype(Psph)
    elseif Psph == Incident() && Pdip == Incident()
        return Incident()
    end
    throw(ErrorException("Cannot perform conversion with given PropagationTypes."))
end
function changerepresentation(
    Tnew::Type{SphericalWaveExpansion{Psph,H,C}},
    dipoles::DipoleArray{Pdip,E,T,C};
    ϵ = 1e-7,
) where {Psph,C,H,Pdip,E,T}
    # Pdual= _dualtype(P)
    Ptmp = _outputmode_dipo2sph(Psph(), Pdip())
    k0 = getwavenumber(dipoles)
    # rsph = (Psph() == Radiated()) ? (2 * _rmax(dipoles)) : (2 * _rmin(dipoles))
    L = equivalentorder(dipoles; ϵ = ϵ)
    Jmax = sℓm_to_j(2, L, L)
    tempcoeffs = _dipole_spherical_innerprod(dipoles, Jmax, Ptmp, k0)
    coefficients = zeros(C, Jmax)
    for (j, val) in enumerate(tempcoeffs)
        s, ℓ, m = j_to_sℓm(j)
        coefficients[sℓm_to_j(s, ℓ, -m)] = (-1)^(m + 1) * val
    end
    return Tnew(SphericalCoefficients(coefficients), k0)
end
function changerepresentation(
    Tnew::Type{SphericalWaveExpansion{Psph}},
    dipoles::DipoleArray{Pdip,E,T,C};
    ϵ = 1e-7,
) where {Psph,C,Pdip,E,T}
    return changerepresentation(
        SphericalWaveExpansion{Psph,SphericalCoefficients{C},C},
        dipoles,
        ϵ = ϵ,
    )
end

function _dipole_spherical_innerprod(
    dipoles::HertzArray{T,C},
    Jmax::Integer,
    P::PropagationType,
    k0::Number,
) where {C,T}
    tempcoeffs = zeros(C, Jmax)
    for (i, position) in enumerate(dipoles.positions)
        Fx, Fy, Fz = F_sℓm_cartesian_array(Jmax, P, position, k0)
        tempcoeffs .+=
            [Fx Fy Fz] * dipoles.orientations[i] * dipoles.dipolemoments[i] * k0 * sqrt(Z₀)
    end
    return tempcoeffs
end
function _dipole_spherical_innerprod(
    dipoles::FitzgeraldArray{T,C},
    Jmax::Integer,
    P::PropagationType,
    k0::Number,
) where {C,T}
    tempcoeffs = zeros(C, Jmax)
    for (i, position) in enumerate(dipoles.positions)
        Fx, Fy, Fz = curlF_sℓm_cartesian_array(Jmax, P, position, k0)
        tempcoeffs .+=
            [Fx Fy Fz] *
            dipoles.orientations[i] *
            dipoles.dipolemoments[i] *
            complex(0.0, -k0) / sqrt(Z₀)
    end
    return tempcoeffs
end

"""
αinc_planewave(L::Integer)


Return incident spherical mode coefficients up to mode order L
of an x-polarized plane wave traveling in negative z-direction.
Compare Hansen: "Spherical Near-Field Measurements" Appendix A1.6 
"""
function αinc_planewave(L::Integer)
    αin = FirstOrderSphericalCoefficients(zeros(ComplexF64, sℓm_to_j(2, L, L)))
    for ℓ = 1:L

        sqrtfac = sqrt((2 * ℓ + 1) / pi)

        αin[sℓm_to_j(1, ℓ, 1)] = -(1im)^(ℓ) * sqrtfac
        αin[sℓm_to_j(2, ℓ, 1)] = -(1im)^(ℓ) * sqrtfac
        αin[sℓm_to_j(1, ℓ, -1)] = -(1im)^(ℓ) * sqrtfac
        αin[sℓm_to_j(2, ℓ, -1)] = (1im)^(ℓ) * sqrtfac
    end
    αin .*= sqrt(Z₀) * 0.5
    return αin
end

"""
αinc_dipole(z::Real, L::Integer, , wavenumber::Real)


Return incident spherical mode coefficients up to mode order L
of an x-polarized Hertzian dipole at `z`.
Compare Hansen: "Spherical Near-Field Measurements" Appendix A1.6 
"""
function αinc_dipole(z::Real, L::Integer, wavenumber::Real)
    αin = FirstOrderSphericalCoefficients(zeros(ComplexF64, sℓm_to_j(2, L, L)))

    hℓ, dhℓ = R_dependencies_array(Radiated(), L, wavenumber * z)
    for ℓ = 1:L
        sqrtfac = sqrt((2 * ℓ + 1) / pi)

        αin[sℓm_to_j(1, ℓ, 1)] = -(1im) * sqrtfac * hℓ[ℓ]
        αin[sℓm_to_j(2, ℓ, 1)] = -(1im) * sqrtfac * hℓ[ℓ]
        αin[sℓm_to_j(1, ℓ, -1)] = sqrtfac * dhℓ[ℓ]
        αin[sℓm_to_j(2, ℓ, -1)] = -sqrtfac * dhℓ[ℓ]
    end
    αin .*= sqrt(Z₀) * 0.25 * wavenumber
    return αin
end

function equivalentorder(coefficients::AbstractSphericalCoefficients; ϵ = 1e-7)
    _, L, __ = j_to_sℓm(length(coefficients))
    return L
end
function equivalentorder(swe::SphericalWaveExpansion; ϵ = 1e-7)
    return equivalentorder(swe.coefficients, ϵ = ϵ)
end

