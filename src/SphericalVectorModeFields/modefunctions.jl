# const Z₀ = 376.730313669

# TODO: Separate spherical coefficients from spherical expansions () the latter are an antenna representation
"""
SphericalExpansions can be Radiating, Incident, or Absorbed
"""
abstract type AbstractSphericalExpansion <: AntennaFieldRepresentation end

function Base.getindex(α::AbstractSphericalExpansion,s::Integer,ℓ::Integer,m::Integer)
    s<1 && DomainError("s must be 1 or 2.")
    s>2 && DomainError("s must be 1 or 2.")
    ℓ<1 && DomainError("ℓ must be positive.")
    abs(m) > ℓ && DomainError("|m| must be ≤ ℓ.")

    return α.coefficients[sℓm_to_j(s,ℓ,m)]
end
function Base.getindex(α::AbstractSphericalExpansion,j::Integer)
    return α.coefficients[j]
end

struct FirstOrder{T<: AbstractSphericalExpansion, C<: Complex}
    coeffs_s1_m_plus::Vector{C}
    coeffs_s1_m_minus::Vector{C}
    coeffs_s2_m_plus::Vector{C}
    coeffs_s2_m_minus::Vector{C}
end

function FirstOrder{T}(α::T) where{T<: AbstractSphericalExpansion}
    _,L,__ = j_to_sℓm(length(α.coefficients))
    coeffs_s1_m_plus = α[1, 1:L, 1]
    coeffs_s1_m_minus = α[1, 1:L, -1]
    coeffs_s2_m_plus = α[2, 1:L, 1]
    coeffs_s2_m_minus = α[2, 1:L, -1]  
    return FirstOrder{T}(coeffs_s1_m_plus, coeffs_s1_m_minus, coeffs_s2_m_plus, coeffs_s2_m_minus)
end
function Base.getindex(α::FirstOrder{T,C}, s::Integer,ℓ::Integer,m::Integer) where{T<: AbstractSphericalExpansion, C<:Complex}
    s<1 && DomainError("s must be 1 or 2.")
    s>2 && DomainError("s must be 1 or 2.")
    ℓ<1 && DomainError("ℓ must be positive.")
    abs(m) > ℓ && DomainError("|m| must be ≤ ℓ.")

    abs(m) != 1 && return C(0)
    if m== 1 
        if s==1 
            return α.coeffs_s1_m_plus
        else return α.coeffs_s2_m_plus
        end
    else 
        if s==1 
            return α.coeffs_s1_m_minus
        else return α.coeffs_s2_m_minus
        end
    end

end
function Base.getindex(α::FirstOrder{T,C}, j::Integer) where{T<: AbstractSphericalExpansion, C<:Complex}
    return α[j_to_sℓm(j)]
end
function Base.length(α::FirstOrder{T,C})  where{T<: AbstractSphericalExpansion, C<:Complex}
    L=length(α.coeffs_s1_m_plus)
    return sℓm_to_j(2,L,L)
end


"""
    RadiatingSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
Store spherical vector wave expansion coefficients for an expansion of radiating type

# Fields 
- `coefficients::Array{C,1}` : array storing the expansion coefficients
"""
mutable struct RadiatingSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
    coefficients::Array{C,1}
end

function RadiatingSphericalExpansion{C}(α::FirstOrder{RadiatingSphericalExpansion{C}}) where C
    L= length(α.coeffs_s1_m_minus)
    J= sℓm_to_j(2,L,L)
    coefficients=zeros(eltype(α.coeffs_s1_m_minus), J)
    for ℓ = 1:L
        coefficients[sℓm_to_j(1, ℓ, 1)] = α.coeffs_s1_m_plus
        coefficients[sℓm_to_j(2, ℓ, 1)] = α.coeffs_s2_m_plus
        coefficients[sℓm_to_j(1, ℓ, -1)] = α.coeffs_s1_m_minus
        coefficients[sℓm_to_j(2, ℓ, -1)] = α.coeffs_s2_m_minus
    end
    return RadiatingSphericalExpansion{C}(coefficients)
end



"""
    IncidentSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
Store spherical vector wave expansion coefficients for an expansion of incident type

# Fields 
- `coefficients::Array{C,1}` : array storing the expansion coefficients
"""
mutable struct IncidentSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
    coefficients::Array{C,1}
end

function IncidentSphericalExpansion{C}(α::FirstOrder{IncidentSphericalExpansion{C}}) where C
    L= length(α.coeffs_s1_m_minus)
    J= sℓm_to_j(2,L,L)
    coefficients=zeros(eltype(α.coeffs_s1_m_minus), J)
    for ℓ = 1:L
        coefficients[sℓm_to_j(1, ℓ, 1)] = α.coeffs_s1_m_plus
        coefficients[sℓm_to_j(2, ℓ, 1)] = α.coeffs_s2_m_plus
        coefficients[sℓm_to_j(1, ℓ, -1)] = α.coeffs_s1_m_minus
        coefficients[sℓm_to_j(2, ℓ, -1)] = α.coeffs_s2_m_minus
    end
    return IncidentSphericalExpansion{C}(coefficients)
end

"""
    AbsorbedSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
Store spherical vector wave expansion coefficients for an expansion of absorbed type

# Fields 
- `coefficients::Array{C,1}` : array storing the expansion coefficients
"""
mutable struct AbsorbedSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
    coefficients::Array{C,1}
end

function AbsorbedSphericalExpansion{C}(α::FirstOrder{AbsorbedSphericalExpansion{C}}) where C
    L= length(α.coeffs_s1_m_minus)
    J= sℓm_to_j(2,L,L)
    coefficients=zeros(eltype(α.coeffs_s1_m_minus), J)
    for ℓ = 1:L
        coefficients[sℓm_to_j(1, ℓ, 1)] = α.coeffs_s1_m_plus
        coefficients[sℓm_to_j(2, ℓ, 1)] = α.coeffs_s2_m_plus
        coefficients[sℓm_to_j(1, ℓ, -1)] = α.coeffs_s1_m_minus
        coefficients[sℓm_to_j(2, ℓ, -1)] = α.coeffs_s2_m_minus
    end
    return AbsorbedSphericalExpansion{C}(coefficients)
end

"""
    UnorthodoxSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
Store spherical vector wave expansion coefficients for an expansion of absorbed type

# Fields 
- `coefficients::Array{C,1}` : array storing the expansion coefficients
"""
mutable struct UnorthodoxSphericalExpansion{C<:Complex} <: AbstractSphericalExpansion
    coefficients::Array{C,1}
end

function UnorthodoxSphericalExpansion{C}(α::FirstOrder{UnorthodoxSphericalExpansion{C}}) where C
    L= length(α.coeffs_s1_m_minus)
    J= sℓm_to_j(2,L,L)
    coefficients=zeros(eltype(α.coeffs_s1_m_minus), J)
    for ℓ = 1:L
        coefficients[sℓm_to_j(1, ℓ, 1)] = α.coeffs_s1_m_plus
        coefficients[sℓm_to_j(2, ℓ, 1)] = α.coeffs_s2_m_plus
        coefficients[sℓm_to_j(1, ℓ, -1)] = α.coeffs_s1_m_minus
        coefficients[sℓm_to_j(2, ℓ, -1)] = α.coeffs_s2_m_minus
    end
    return UnorthodoxSphericalExpansion{C}(coefficients)
end

function converttype(T::Type{RadiatingSphericalExpansion{C}}, α::AbstractSphericalExpansion) where {C<:Complex}
    return RadiatingSphericalExpansion(α.coefficients)
end
function converttype(T::Type{IncidentSphericalExpansion{C}}, α::AbstractSphericalExpansion) where {C<:Complex}
    return IncidentSphericalExpansion(α.coefficients)
end
function converttype(T::Type{AbsorbedSphericalExpansion{C}}, α::AbstractSphericalExpansion) where {C<:Complex}
    return AbsorbedSphericalExpansion(α.coefficients)
end
function converttype(T::Type{UnorthodoxSphericalExpansion{C}}, α::AbstractSphericalExpansion) where {C<:Complex}
    return UnorthodoxSphericalExpansion(α.coefficients)
end

function elementtype(T::Type{<:AbstractSphericalExpansion})
    return ComplexF64
end
function elementtype(T::Type{RadiatingSphericalExpansion{C}}) where {C<:Complex}
    return C
end
function elementtype(T::Type{IncidentSphericalExpansion{C}}) where {C<:Complex}
    return C
end
function elementtype(T::Type{AbsorbedSphericalExpansion{C}}) where {C<:Complex}
    return C
end
function elementtype(T::Type{UnorthodoxSphericalExpansion{C}}) where {C<:Complex}
    return C
end

function elementtype(a::AbstractSphericalExpansion)
    return ComplexF64
end
function elementtype(a::RadiatingSphericalExpansion{C}) where {C<:Complex}
    return C
end
function elementtype(a::IncidentSphericalExpansion{C}) where {C<:Complex}
    return C
end
function elementtype(a::AbsorbedSphericalExpansion{C}) where {C<:Complex}
    return C
end
function elementtype(a::UnorthodoxSphericalExpansion{C}) where {C<:Complex}
    return C
end


function reciprocaltype(T::Type{RadiatingSphericalExpansion{C}}) where {C<:Complex}
    return IncidentSphericalExpansion{C}
end
function reciprocaltype(T::Type{IncidentSphericalExpansion{C}}) where {C<:Complex}
    return RadiatingSphericalExpansion{C}
end
# function reciprocaltype(a::RadiatingSphericalExpansion{C}) where{C<:Complex}
#     return IncidentSphericalExpansion{C}
# end
# function reciprocaltype(a::IncidentSphericalExpansion{C}) where{C<:Complex}
#     return RadiatingSphericalExpansion{C}
# end

"""
    sℓm_to_j(s,ℓ,m)

    Convert multi-index s ℓ m to single index j
"""
function sℓm_to_j(s, ℓ, m)
    j = 2 * (ℓ * (ℓ + 1) + m - 1) + s
    return j

end

"""
    j_to_sℓm(j)

    Convert single index j to multi-index s ℓ m
"""
function j_to_sℓm(j::Integer)
    jtype=typeof(j)
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
    F_sℓm_spherical_array(Jmaxxc,r,ϑ,φ, k0) -> Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
  
Return all spherical vector wave functions upt to Jmaxx at r, ϑ, φ in spherical coordinates
"""
function F_sℓm_spherical_array(Jmaxx::Integer, T::Type{<:AbstractSphericalExpansion}, r::Real, ϑ::Number, φ::Real, k0::Real)
    _, Lmax, __ = j_to_sℓm(Jmaxx)
    Jmax = 2 * Lmax * (Lmax + 2)
    kA = k0 * r

    Fr = zeros(elementtype(T), Jmax)
    Fϑ = zeros(elementtype(T), Jmax)
    Fφ = zeros(elementtype(T), Jmax)

    if abs(kA < 100 * eps())
        j = sℓm_to_j(2, 1, -1)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, -1, T, ϑ, φ, k0)
        j = sℓm_to_j(2, 1, 0)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, 0, T, ϑ, φ, k0)
        j = sℓm_to_j(2, 1, 1)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, 1, T, ϑ, φ, k0)

        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    end

    if abs(ϑ) <= 100 * eps()
        for ℓ in 1:Lmax
            j = sℓm_to_j(1, ℓ, -1)
            Fr[j], Fϑ[j], Fφ[j] = F_sℓm_thetazero(1, ℓ, -1, T, r, φ, k0)
            Fr[j + 1], Fϑ[j + 1], Fφ[j + 1] = F_sℓm_thetazero(2, ℓ, -1, T, r, φ, k0)
            Fr[j + 2], Fϑ[j + 2], Fφ[j + 2] = F_sℓm_thetazero(1, ℓ, 0, T, r, φ, k0)
            Fr[j + 3], Fϑ[j + 3], Fφ[j + 3] = F_sℓm_thetazero(2, ℓ, 0, T, r, φ, k0)
            Fr[j + 4], Fϑ[j + 4], Fφ[j + 4] = F_sℓm_thetazero(1, ℓ, 1, T, r, φ, k0)
            Fr[j + 5], Fϑ[j + 5], Fφ[j + 5] = F_sℓm_thetazero(2, ℓ, 1, T, r, φ, k0)
        end
        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    elseif (abs(ϑ - pi) <= 100 * eps())
        for ℓ in 1:Lmax
            j = sℓm_to_j(1, ℓ, -1)
            Fr[j], Fϑ[j], Fφ[j] = F_sℓm_thetapi(1, ℓ, -1, T, r, φ, k0)
            Fr[j + 1], Fϑ[j + 1], Fφ[j + 1] = F_sℓm_thetapi(2, ℓ, -1, T, r, φ, k0)
            Fr[j + 2], Fϑ[j + 2], Fφ[j + 2] = F_sℓm_thetapi(1, ℓ, 0, T, r, φ, k0)
            Fr[j + 3], Fϑ[j + 3], Fφ[j + 3] = F_sℓm_thetapi(2, ℓ, 0, T, r, φ, k0)
            Fr[j + 4], Fϑ[j + 4], Fφ[j + 4] = F_sℓm_thetapi(1, ℓ, 1, T, r, φ, k0)
            Fr[j + 5], Fϑ[j + 5], Fφ[j + 5] = F_sℓm_thetapi(2, ℓ, 1, T, r, φ, k0)
        end
        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    end

    zℓ = zeros(elementtype(T), Lmax)
    dzℓ = zeros(elementtype(T), Lmax)
    for ℓ in eachindex(zℓ)
        oneoversqrt = 1 / (sqrt(ℓ * (ℓ + 1)))
        zℓ[ℓ] = zc_ℓ(T, ℓ, kA) * oneoversqrt
        dzℓ[ℓ] = oneoverkA_deriv_zc_ℓ(T, ℓ, kA) * oneoversqrt
    end


    for m in 1:Lmax
        P1, P2, P3 = legendre_deps_array(m, Lmax, ϑ)
        exp_jmφ = cis(m * φ)  # exp(1im*m*φ)
        exp_minusjmφ = conj(exp_jmφ) * (-1)^m
        for ℓ in m:Lmax
            ℓindex = ℓ - m + 1
            j = sℓm_to_j(1, ℓ, -m)
            Fϑ[j] = -1im * zℓ[ℓ] * P1[ℓindex] * exp_minusjmφ
            Fφ[j] = -zℓ[ℓ] * P2[ℓindex] * exp_minusjmφ

            j = sℓm_to_j(2, ℓ, -m)
            Fr[j] = ℓ * (ℓ + 1) / kA * zℓ[ℓ] * P3[ℓindex] * exp_minusjmφ
            Fϑ[j] = dzℓ[ℓ] * P2[ℓindex] * exp_minusjmφ
            Fφ[j] = -1im * dzℓ[ℓ] * P1[ℓindex] * exp_minusjmφ

            j = sℓm_to_j(1, ℓ, m)
            Fϑ[j] = 1im * zℓ[ℓ] * P1[ℓindex] * exp_jmφ
            Fφ[j] = -zℓ[ℓ] * P2[ℓindex] * exp_jmφ

            j = sℓm_to_j(2, ℓ, m)
            Fr[j] = ℓ * (ℓ + 1) / kA * zℓ[ℓ] * P3[ℓindex] * exp_jmφ
            Fϑ[j] = dzℓ[ℓ] * P2[ℓindex] * exp_jmφ
            Fφ[j] = 1im * dzℓ[ℓ] * P1[ℓindex] * exp_jmφ
        end
    end
    P1, P2, P3 = legendre_deps_array(0, Lmax, ϑ)
    for ℓ in 1:Lmax # m==0
        j = sℓm_to_j(1, ℓ, 0)
        Fϑ[j] = 1im * zℓ[ℓ] * P1[ℓ + 1]
        Fφ[j] = -zℓ[ℓ] * P2[ℓ + 1]

        j = sℓm_to_j(2, ℓ, 0)
        Fr[j] = ℓ * (ℓ + 1) / kA * zℓ[ℓ] * P3[ℓ + 1]
        Fϑ[j] = dzℓ[ℓ] * P2[ℓ + 1]
        Fφ[j] = 1im * dzℓ[ℓ] * P1[ℓ + 1]

    end

    return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
end

"""
    F_sℓm_thetazero(s,ℓ,m,c,r,φ,k0)- > Fr,Fϑ,Fφ

Catches special cases for ϑ==0
Hansen p. 3245f.
"""
function F_sℓm_thetazero(s::Integer, ℓ::Integer, m::Integer, T::Type{<:AbstractSphericalExpansion}, r::Real, φ::Real, k0::Real)
    Fr = zero(elementtype(T))
    Fϑ = zero(elementtype(T))
    Fφ = zero(elementtype(T))
    if abs(m) > 1
        return [0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im]

    elseif s == 1
        fac = -sqrt((2 * ℓ + 1) / pi) / 4 * zc_ℓ(T, ℓ, k0 * r) * 1im * cis(φ * m)
        Fϑ = abs(m) * fac
        Fφ = 1im * m * fac
    elseif s == 2
        if m == 0
            Fr = sqrt(ℓ * (ℓ + 1) * (2 * ℓ + 1) / (4 * pi)) * zc_ℓ(T, ℓ, k0 * r) / (k0 * r)
        else
            fac = -m * sqrt((2 * ℓ + 1) / pi) / 4 * oneoverkA_deriv_zc_ℓ(T, ℓ, k0 * r) * cis(m * φ)
            Fϑ = fac
            Fφ = 1im * m * fac
        end
    end
    return Fr, Fϑ, Fφ
end

"""
    F_sℓm_thetapi(s,ℓ,m,c,r,φ,k0) -> Fr,Fϑ,Fφ

Catches special cases for ϑ==pi
Hansen p. 3245f.
"""
function F_sℓm_thetapi(s::Integer, ℓ::Integer, m::Integer, T::Type{<:AbstractSphericalExpansion}, r::Real, φ::Real, k0::Real)
    Fr, Fϑ, Fφ = F_sℓm_thetazero(s, ℓ, m, T, r, φ, k0)
    Fr *= (-1)^ℓ
    Fϑ *= (-1)^(ℓ + s)
    Fφ *= (-1)^(ℓ + s + 1)
    return Fr, Fϑ, Fφ
end

"""
    F_sℓm_rzero(s,ℓ, m ,c,ϑ,φ, k0) -> Fr,Fϑ,Fφ

Catches special cases for r==0
Hansen p. 3245f.
"""
function F_sℓm_spherical_rzero(s::Integer, ℓ::Integer, m::Integer, T::Type{<:AbstractSphericalExpansion}, ϑ::Number, φ::Real, k0::Real)
    return convert(elementtype(T), Inf), convert(elementtype(T), Inf), convert(elementtype(T), Inf)
end
function F_sℓm_spherical_rzero(
    s::Integer, ℓ::Integer, m::Integer, T::Type{IncidentSphericalExpansion{C}}, ϑ::Number, φ::Real, k0::Real
) where {C<:Complex}

    Fr = zero(elementtype(T))
    Fϑ = zero(elementtype(T))
    Fφ = zero(elementtype(T))
    if (ℓ == 1) && s == 2
        if m == 0
            fac = sqrt(6 / pi) / 6
            Fr = fac * cos(ϑ)
            Fϑ = -fac * sin(ϑ)
        elseif abs(m) == 1
            fac = -m * sqrt(3 / pi) / 6 * cis(m * φ)
            Fr = fac * sin(ϑ)
            Fϑ = fac * cos(ϑ)
            Fφ = complex(0.0, m) * fac
        end
    end
    return Fr, Fϑ, Fφ
end


"""
    F_sℓm_spherical(s,ℓ, m,r,ϑ,φ) -> [Fr;Fϑ;Fφ]

Return normalized vector spherical wave function at r,ϑ,φ in spherical coordinates 
"""
function F_sℓm_spherical(
    s::Integer, ℓ::Integer, m::Integer, T::Type{<:AbstractSphericalExpansion}, r::Real, ϑ::Number, φ::Real, k0::Real
)
    j = sℓm_to_j(s, ℓ, m)
    Fr, Fϑ, F_sℓm_spherical_array(j, T, r, ϑ, φ, k0)
    return [Fr[j]; Fϑ[j]; Fφ[j]]
end

"""
    F_sℓm_cartesian(s,ℓ,m,c,R,k0) -> [Fx;Fy;Fz]

Return normalized vector spherical wave function at R in cartesian coordinates 
"""
function F_sℓm_cartesian(s, ℓ, m, T, R, k0)
    r = norm(R)
    φ = mod(atan(R[2], R[1]), 2pi)
    ϑ = abs(mod(atan(sqrt(R[1]^2 + R[2]^2), R[3]) + pi, 2pi) - pi)
    Fspherical = F_sℓm_spherical(s, ℓ, m, T, r, ϑ, φ, k0)
    sint, cost = sincos(ϑ)
    sinp, cosp = sincos(φ)
    Fcartesian = [
        sint*cosp cost*cosp -sinp
        sint*sinp cost*sinp cosp
        cost -sint 0
    ] * Fspherical
    return Fcartesian
end

"""
    F_sℓm_cartesian_array(s,ℓ,m,c,R,k0) -> [Fx;Fy;Fz]

Return normalized vector spherical wave function at R in cartesian coordinates 
"""
function F_sℓm_cartesian_array(
    Jmaxx::Integer, T::Type{<:AbstractSphericalExpansion}, R::AbstractArray{N}, k0::Number
) where {N<:Number}
    C = elementtype(T)
    r = norm(R)
    φ = zero(N)
    ϑ = zero(N)
    if k0 * r > 100 * eps()
        φ = (mod(atan(R[2], R[1]), 2pi))
        ϑ = abs(mod(atan(sqrt(R[1]^2 + R[2]^2), R[3]) + pi, 2pi) - pi)
        # if ϑ > pi
        #     ϑ=(2*pi-ϑ)
        # end 
    end
    Fr, Fϑ, Fφ = F_sℓm_spherical_array(Jmaxx, T, r, ϑ, φ, k0)
    sint, cost = sincos(ϑ)
    sinp, cosp = sincos(φ)
    Fcartesian = [Fr Fϑ Fφ] * [
                       sint*cosp   sint*sinp    cost
        cost*cosp   cost*sinp   -sint
        -sinp cosp 0
    ]
    return convert.(C, Fcartesian[:, 1]), convert.(C, Fcartesian[:, 2]), convert.(C, Fcartesian[:, 3])
end

function curlF_sℓm_cartesian_array(
    Jmaxx::Integer, T::Type{<:AbstractSphericalExpansion}, R::AbstractArray{N}, k0::Number
) where {N<:Number}
    Fx, Fy, Fz = F_sℓm_cartesian_array(Jmaxx, T, R, k0)

    Fx[1:2:end], Fx[2:2:end] = Fx[2:2:end], Fx[1:2:end]
    Fy[1:2:end], Fy[2:2:end] = Fy[2:2:end], Fy[1:2:end]
    Fz[1:2:end], Fz[2:2:end] = Fz[2:2:end], Fz[1:2:end]

    return Fx, Fy, Fz
end

"""
    K_sℓm_thetazero(s,ℓ,m,c,φ) -> Kϑ, Kφ

Catches special cases for ϑ==0
Hansen p. 329f.
"""
function K_sℓm_thetazero(s::Integer, ℓ::Integer, m::Integer, φ::Real, outputtype::Type{C}=ComplexF64) where {C<:Complex}
    Kϑ = zero(outputtype)
    Kφ = zero(outputtype)

    if abs(m) > 1
        return zero(outputtype), zero(outputtype)

    elseif s == 1
        fac = outputtype(-sqrt((2 * ℓ + 1) / pi) / 4 * (1im)^(ℓ + 1) * 1im * cis(φ * m))
        Kϑ = abs(m) * fac
        Kφ = outputtype(0, 1) * m * fac
        return Kϑ, Kφ
    elseif s == 2
        fac = outputtype(-m * sqrt((2 * ℓ + 1) / pi) / 4 * (1im)^(ℓ) * cis(m * φ))
        Kϑ = fac
        Kφ = outputtype(0, 1) * m * fac
        return Kϑ, Kφ
    end

    return Kϑ, Kφ
end


"""
    K_sℓm_thetapi(s,ℓ,m,c,φ) -> Kϑ, Kφ

Catches special cases for ϑ==π
Hansen p. 329f.
"""
function K_sℓm_thetapi(s::Integer, ℓ::Integer, m::Integer, φ::Real)
    Kϑ, Kφ = K_sℓm_thetazero(s, ℓ, m, φ)
    Kϑ *= (-1)^(ℓ + s)
    Kφ *= (-1)^(ℓ + s - 1)
    return Kϑ, Kφ
end

function ehfield(α::AbstractSphericalExpansion, R::AbstractVector{<:Number}, k0::Number)
    C = elementtype(typeof(α))
    sqrtZ₀ = C(sqrt(Z₀))

    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(α.coefficients)

        Fx, Fy, Fz = F_sℓm_cartesian_array(J, typeof(α), R, k0)
        Ecartesian = SVector{3,C}((C(k0) * sqrtZ₀) .* (udot(F, α.coefficients) for F in (Fx, Fy, Fz)))

        Hx = udot(Fx[1:2:J], α.coefficients[2:2:J]) + udot(Fx[2:2:J], α.coefficients[1:2:J])
        Hy = udot(Fy[1:2:J], α.coefficients[2:2:J]) + udot(Fy[2:2:J], α.coefficients[1:2:J])
        Hz = udot(Fz[1:2:J], α.coefficients[2:2:J]) + udot(Fz[2:2:J], α.coefficients[1:2:J])
        Hcartesian = SVector{3,C}((C(0.0, k0) / sqrtZ₀) .* (Hx, Hy, Hz))

        return Ecartesian, Hcartesian

    else
        Ecartesian = _E_at_origin(α::IncidentSphericalExpansion{C}, k0::Number)::SVector{3,C}
        Hcartesian = _H_at_origin(α::IncidentSphericalExpansion{C}, k0::Number)::SVector{3,C}
        return Ecartesian, Hcartesian
    end
end

function efield(α::AbstractSphericalExpansion, R::AbstractVector{<:Number}, k0::Number)
    C = elementtype(typeof(α))
    sqrtZ₀ = convert(C, sqrt(Z₀))

    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(α.coefficients)

        Fx, Fy, Fz = F_sℓm_cartesian_array(J, typeof(α), R, k0)

        return SVector{3,C}((C(k0) * sqrtZ₀) .* (udot(F, α.coefficients) for F in (Fx, Fy, Fz)))

    else
        return _E_at_origin(α, k0::Number)::SVector{3,C}
    end
end

function hfield(α::AbstractSphericalExpansion, R::AbstractVector{<:Number}, k0::Number)
    C = elementtype(typeof(α))
    sqrtZ₀ = convert(C, sqrt(Z₀))

    ϵ = 1e-10
    if k0 * norm(R) > ϵ
        J = length(α.coefficients)

        Fx, Fy, Fz = curlF_sℓm_cartesian_array(J, typeof(α), R, k0)

        return SVector{3,C}((C(0.0, k0) / sqrtZ₀) .* (udot(F, α.coefficients) for F in (Fx, Fy, Fz)))

    else
        return _H_at_origin(α, k0::Number)::SVector{3,C}
    end
end

function _F2m1cartesian_at_origin(C::Type{<:Complex}=ComplexF64)
    F2m11 = SVector{3}(C(√(3 / pi) / 6), C(-1im * √(3 / pi) / 6), C(0.0))
    F201 = SVector{3}(C(0.0), C(0.0), C(√(6 / pi) / 6))
    F211 = SVector{3}(C(-√(3 / pi) / 6), C(-1im * √(3 / pi) / 6), C(0.0))
    return F2m11, F201, F211
end

function _H_at_origin(α::AbstractSphericalExpansion, k0::Number)
    C = elementtype(typeof(α))
    return SVector{3}(_H_at_origin(converttype(IncidentSphericalExpansion{C}, α), k0) .+ C(Inf, 0.0))
end
function _H_at_origin(α::IncidentSphericalExpansion{C}, k0::Number) where {C<:Complex}
    F2m11, F201, F211 = _F2m1cartesian_at_origin(C)
    Hcartesian = SVector{3}(
        C(0.0, k0) / C(sqrt(Z₀)) *
        C.(
            F211 * α.coefficients[sℓm_to_j(1, 1, 1)] +
            F2m11 * α.coefficients[sℓm_to_j(1, 1, -1)] +
            F201 * α.coefficients[sℓm_to_j(1, 1, 0)]
        ),
    )
    return Hcartesian

end

function _E_at_origin(α::AbstractSphericalExpansion, k0::Number)
    C = elementtype(typeof(α))
    return SVector{3}(_E_at_origin(converttype(IncidentSphericalExpansion{C}, α), k0) .+ C(0.0, Inf))
end
function _E_at_origin(α::IncidentSphericalExpansion{C}, k0::Number) where {C<:Complex}
    F2m11, F201, F211 = _F2m1cartesian_at_origin(C)
    Ecartesian = SVector{3}(
        C(k0) *
        C(sqrt(Z₀)) *
        C.(
            F211 * α.coefficients[sℓm_to_j(2, 1, 1)] +
            F2m11 * α.coefficients[sℓm_to_j(2, 1, -1)] +
            F201 * α.coefficients[sℓm_to_j(2, 1, 0)]
        ),
    )
    return Ecartesian
end

"""
    K_sℓm(s,ℓ,m,ϑ,φ) -> Kϑ,Kφ

Return farfield vector spherical wave function in spherical coordinates 
"""
function K_sℓm(s::Integer, ℓ::Integer, m::Integer, ϑ::Number, φ::Real, outputtype::Type{C}=ComplexF64) where {C<:Complex}
    j = sℓm_to_j(s, ℓ, m)
    Kϑ, Kφ = K_sℓm_array(j, ϑ, φ, outputtype)
    return Kϑ[end], Kφ[end]
end

"""
    K_sℓm_array(Jmaxx,ϑ,φ) -> Kϑ[1:Jmaxx], Kφ[1:Jmaxx]

Return farfield for all vector spherical wave functions up to Jmaxx in spherical coordinates
"""
function K_sℓm_array(Jmaxx::Integer, ϑ::Number, φ::Real, outputtype::Type{C}=ComplexF64) where {C<:Complex}
    Kϑ, Kφ = K_sℓm_array(Jmaxx::Integer, ϑ::Number, [φ], outputtype)
    return Kϑ[:, 1], Kφ[:, 1]
end
function K_sℓm_array(Jmaxx::Integer, ϑ::Number, φvec::Array{<:Real,1}, outputtype::Type{C}=ComplexF64) where {C<:Complex}
    _, Lmax, __ = j_to_sℓm(Jmaxx)
    Jmax = 2 * Lmax * (Lmax + 2)
    Nφ = length(φvec)

    Kϑ = zeros(outputtype, Jmax, Nφ)
    Kφ = zeros(outputtype, Jmax, Nφ)
    if (abs(ϑ) > 10 * eps()) && (abs(ϑ - pi) > 10 * eps())
        zℓ = [outputtype((1im)^(ℓ + 1) / (sqrt(ℓ * (ℓ + 1)))) for ℓ in 1:Lmax]
        dzℓ = [outputtype((1im)^(ℓ) / (sqrt(ℓ * (ℓ + 1)))) for ℓ in 1:Lmax]

        for m in 1:Lmax
            P1, P2, __ = legendre_deps_array(m, Lmax, ϑ)
            exp_jmφ = outputtype.(cis.(m * φvec))
            exp_minusjmφ = conj(exp_jmφ) * (-1)^m
            for ℓ in m:Lmax
                ℓindex = ℓ - m + 1
                j = sℓm_to_j(1, ℓ, -m)
                Kϑ[j, :] = -1im * zℓ[ℓ] * P1[ℓindex] * exp_minusjmφ
                Kφ[j, :] = -zℓ[ℓ] * P2[ℓindex] * exp_minusjmφ

                j = sℓm_to_j(2, ℓ, -m)
                Kϑ[j, :] = dzℓ[ℓ] * P2[ℓindex] * exp_minusjmφ
                Kφ[j, :] = -1im * dzℓ[ℓ] * P1[ℓindex] * exp_minusjmφ

                j = sℓm_to_j(1, ℓ, m)
                Kϑ[j, :] = 1im * zℓ[ℓ] * P1[ℓindex] * exp_jmφ
                Kφ[j, :] = -zℓ[ℓ] * P2[ℓindex] * exp_jmφ

                j = sℓm_to_j(2, ℓ, m)
                Kϑ[j, :] = dzℓ[ℓ] * P2[ℓindex] * exp_jmφ
                Kφ[j, :] = 1im * dzℓ[ℓ] * P1[ℓindex] * exp_jmφ
            end
        end
        P1, P2, __ = legendre_deps_array(0, Lmax, ϑ)
        for ℓ in 1:Lmax # m==0
            j = sℓm_to_j(1, ℓ, 0)
            Kϑ[j, :] = 1im * zℓ[ℓ] * P1[ℓ + 1] * ones(outputtype, Nφ)
            Kφ[j, :] = -zℓ[ℓ] * P2[ℓ + 1] * ones(outputtype, Nφ)

            j = sℓm_to_j(2, ℓ, 0)
            Kϑ[j, :] = dzℓ[ℓ] * P2[ℓ + 1] * ones(outputtype, Nφ)
            Kφ[j, :] = 1im * dzℓ[ℓ] * P1[ℓ + 1] * ones(outputtype, Nφ)
        end

        return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
    else
        if abs(ϑ) <= 10 * eps()
            for ℓ in 1:Lmax
                j = sℓm_to_j(1, ℓ, -1)
                for (kφ, φ) in enumerate(φvec)
                    Kϑ[j, kφ], Kφ[j, kφ] = K_sℓm_thetazero(1, ℓ, -1, φ)
                    Kϑ[j + 1, kφ], Kφ[j + 1, kφ] = K_sℓm_thetazero(2, ℓ, -1, φ)

                    Kϑ[j + 4, kφ], Kφ[j + 4, kφ] = K_sℓm_thetazero(1, ℓ, 1, φ)
                    Kϑ[j + 5, kφ], Kφ[j + 5, kφ] = K_sℓm_thetazero(2, ℓ, 1, φ)
                end
            end
            return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
        elseif (abs(ϑ - pi) <= 10 * eps())
            for ℓ in 1:Lmax
                j = sℓm_to_j(1, ℓ, -1)
                for (kφ, φ) in enumerate(φvec)
                    Kϑ[j, kφ], Kφ[j, kφ] = K_sℓm_thetapi(1, ℓ, -1, φ)
                    Kϑ[j + 1, kφ], Kφ[j + 1, kφ] = K_sℓm_thetapi(2, ℓ, -1, φ)

                    Kϑ[j + 4, kφ], Kφ[j + 4, kφ] = K_sℓm_thetapi(1, ℓ, 1, φ)
                    Kϑ[j + 5, kφ], Kφ[j + 5, kφ] = K_sℓm_thetapi(2, ℓ, 1, φ)
                end
            end
            return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
        end
        return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
    end


end

function farfield(α::AbstractSphericalExpansion, ϑ, φ)
    return Sphfarfield(α.coefficients, ϑ, φ)
end

function Sphfarfield(α::AbstractVector, ϑ::Number, φ::Real)
    J = length(α)
    Kϑ, Kφ = K_sℓm_array(J, ϑ, φ)
    Fϑ = udot(Kϑ[1:J], α)
    Fφ = udot(Kφ[1:J], α)

    return sqrt(Z₀) * Fϑ, sqrt(Z₀) * Fφ
end

function Sphfarfield(α::AbstractVector, ϑ::Number, φ::Array{<:Real,1})
    C = eltype(α)
    J = length(α)
    Kϑ, Kφ = K_sℓm_array(J, ϑ, φ, C)
    Fϑ = zeros(C, length(φ))
    Fφ = zeros(C, length(φ))
    for k in eachindex(φ)
        Fϑ[k] = udot(Kϑ[1:J, k], α)
        Fφ[k] = udot(Kφ[1:J, k], α)
    end

    return C(sqrt(Z₀)) .* Fϑ, C(sqrt(Z₀)) .* Fφ
end

function Sphfarfield(α::AbstractVector, ϑ::Array{<:Number,1}, φ::Array{<:Real,1})
    C = eltype(α)
    J = length(α)
    Fϑ = zeros(C, length(ϑ), length(φ))
    Fφ = zeros(C, length(ϑ), length(φ))
    for kk in 1:length(ϑ)
        Kϑ, Kφ = K_sℓm_array(J, ϑ[kk], φ, C)
        for k in eachindex(φ)
            Fϑ[kk, k] = udot(Kϑ[1:J, k], α)
            Fφ[kk, k] = udot(Kφ[1:J, k], α)
        end
    end

    return C(sqrt(Z₀)) .* Fϑ, C(sqrt(Z₀)) .* Fφ
end

function transmission(αrad::RadiatingSphericalExpansion, αinc::IncidentSphericalExpansion, ::Number)
    return transmission(αrad, αinc)
end
function transmission(αrad::RadiatingSphericalExpansion, αinc::IncidentSphericalExpansion)
    T = promote_type(elementtype(αrad), elementtype(αinc))
    _, Lrad, __ = j_to_sℓm(length(αrad.coefficients))
    _, Linc, __ = j_to_sℓm(length(αinc.coefficients))
    L = minimum([Linc, Lrad])
    J = sℓm_to_j(2, L, L)

    return -udot(αtoβ(αrad.coefficients[1:J]), αinc.coefficients[1:J])
end

function transmission(αinc::IncidentSphericalExpansion, αrad::RadiatingSphericalExpansion)
    return transmission(αrad, αinc)
end
function transmission(αinc::IncidentSphericalExpansion, αrad::RadiatingSphericalExpansion, ::Number)
    return transmission(αinc, αrad)
end
