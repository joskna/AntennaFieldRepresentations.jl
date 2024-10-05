"""
    F_sℓm_spherical_array(Jmaxxc,r,ϑ,φ, k0) -> Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
  
Return all spherical vector wave functions upt to Jmaxx at r, ϑ, φ in spherical coordinates
"""
function F_sℓm_spherical_array(
    Jmaxx::Integer,
    P:: PropagationType,
    r::T,
    ϑ::Number,
    φ::T,
    k0::T,
) where {T<:Real}
    _, Lmax, __ = j_to_sℓm(Jmaxx)
    Jmax = 2 * Lmax * (Lmax + 2)
    kA = k0 * r

    Fr = zeros(Complex{T}, Jmax)
    Fϑ = zeros(Complex{T}, Jmax)
    Fφ = zeros(Complex{T}, Jmax)

    if abs(kA < 100 * eps())
        j = sℓm_to_j(2, 1, -1)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, -1, P, ϑ, φ, k0)
        j = sℓm_to_j(2, 1, 0)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, 0, P, ϑ, φ, k0)
        j = sℓm_to_j(2, 1, 1)
        Fr[j], Fϑ[j], Fφ[j] = F_sℓm_spherical_rzero(2, 1, 1, P, ϑ, φ, k0)

        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    end

    if abs(ϑ) <= 100 * eps()
        for ℓ = 1:Lmax
            j = sℓm_to_j(1, ℓ, -1)
            Fr[j], Fϑ[j], Fφ[j] = F_sℓm_thetazero(1, ℓ, -1, P, r, φ, k0)
            Fr[j+1], Fϑ[j+1], Fφ[j+1] = F_sℓm_thetazero(2, ℓ, -1, P, r, φ, k0)
            Fr[j+2], Fϑ[j+2], Fφ[j+2] = F_sℓm_thetazero(1, ℓ, 0, P, r, φ, k0)
            Fr[j+3], Fϑ[j+3], Fφ[j+3] = F_sℓm_thetazero(2, ℓ, 0, P, r, φ, k0)
            Fr[j+4], Fϑ[j+4], Fφ[j+4] = F_sℓm_thetazero(1, ℓ, 1, P, r, φ, k0)
            Fr[j+5], Fϑ[j+5], Fφ[j+5] = F_sℓm_thetazero(2, ℓ, 1, P, r, φ, k0)
        end
        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    elseif (abs(ϑ - pi) <= 100 * eps())
        for ℓ = 1:Lmax
            j = sℓm_to_j(1, ℓ, -1)
            Fr[j], Fϑ[j], Fφ[j] = F_sℓm_thetapi(1, ℓ, -1, P, r, φ, k0)
            Fr[j+1], Fϑ[j+1], Fφ[j+1] = F_sℓm_thetapi(2, ℓ, -1, P, r, φ, k0)
            Fr[j+2], Fϑ[j+2], Fφ[j+2] = F_sℓm_thetapi(1, ℓ, 0, P, r, φ, k0)
            Fr[j+3], Fϑ[j+3], Fφ[j+3] = F_sℓm_thetapi(2, ℓ, 0, P, r, φ, k0)
            Fr[j+4], Fϑ[j+4], Fφ[j+4] = F_sℓm_thetapi(1, ℓ, 1, P, r, φ, k0)
            Fr[j+5], Fϑ[j+5], Fφ[j+5] = F_sℓm_thetapi(2, ℓ, 1, P, r, φ, k0)
        end
        return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
    end

    zℓ = zeros(Complex{T}, Lmax)
    dzℓ = zeros(Complex{T}, Lmax)
    for ℓ in eachindex(zℓ)
        oneoversqrt = 1 / (sqrt(ℓ * (ℓ + 1)))
        zℓ[ℓ] = zc_ℓ(P, ℓ, kA) * oneoversqrt
        dzℓ[ℓ] = oneoverkA_deriv_zc_ℓ(P, ℓ, kA) * oneoversqrt
    end


    for m = 1:Lmax
        P1, P2, P3 = legendre_deps_array(m, Lmax, ϑ)
        exp_jmφ = cis(m * φ)  # exp(1im*m*φ)
        exp_minusjmφ = conj(exp_jmφ) * (-1)^m
        for ℓ = m:Lmax
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
    for ℓ = 1:Lmax # m==0
        j = sℓm_to_j(1, ℓ, 0)
        Fϑ[j] = 1im * zℓ[ℓ] * P1[ℓ+1]
        Fφ[j] = -zℓ[ℓ] * P2[ℓ+1]

        j = sℓm_to_j(2, ℓ, 0)
        Fr[j] = ℓ * (ℓ + 1) / kA * zℓ[ℓ] * P3[ℓ+1]
        Fϑ[j] = dzℓ[ℓ] * P2[ℓ+1]
        Fφ[j] = 1im * dzℓ[ℓ] * P1[ℓ+1]

    end

    return Fr[1:Jmaxx], Fϑ[1:Jmaxx], Fφ[1:Jmaxx]
end

"""
    F_sℓm_thetazero(s,ℓ,m,c,r,φ,k0)- > Fr,Fϑ,Fφ

Catches special cases for ϑ==0
Hansen p. 3245f.
"""
function F_sℓm_thetazero(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    P::PropagationType,
    r::T,
    φ::T,
    k0::T,
) where {T<:Real}
    Fr = zero(Complex{T})
    Fϑ = zero(Complex{T})
    Fφ = zero(Complex{T})
    if abs(m) > 1
        return [0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im]

    elseif s == 1
        fac = -sqrt((2 * ℓ + 1) / pi) / 4 * zc_ℓ(P, ℓ, k0 * r) * 1im * cis(φ * m)
        Fϑ = abs(m) * fac
        Fφ = 1im * m * fac
    elseif s == 2
        if m == 0
            Fr = sqrt(ℓ * (ℓ + 1) * (2 * ℓ + 1) / (4 * pi)) * zc_ℓ(P, ℓ, k0 * r) / (k0 * r)
        else
            fac =
                -m * sqrt((2 * ℓ + 1) / pi) / 4 *
                oneoverkA_deriv_zc_ℓ(P, ℓ, k0 * r) *
                cis(m * φ)
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
function F_sℓm_thetapi(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    P::PropagationType,
    r::T,
    φ::T,
    k0::T,
) where{T<: Real}
    Fr, Fϑ, Fφ = F_sℓm_thetazero(s, ℓ, m, P, r, φ, k0)
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
function F_sℓm_spherical_rzero(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    P::PropagationType,
    ϑ::Number,
    φ::T,
    k0::T,
) where{T<: Real}
    return convert(Complex{T}, Inf),
    convert(Complex{T}, Inf),
    convert(Complex{T}, Inf)
end
function F_sℓm_spherical_rzero(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    P::Incident,
    ϑ::Number,
    φ::T,
    k0::T,
) where{T<: Real}

    Fr = zero(Complex{T})
    Fϑ = zero(Complex{T})
    Fφ = zero(Complex{T})
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
    s::Integer,
    ℓ::Integer,
    m::Integer,
    P::PropagationType,
    r::T,
    ϑ::Number,
    φ::T,
    k0::T,
) where{T<: Real}
    j = sℓm_to_j(s, ℓ, m)
    Fr, Fϑ, F_sℓm_spherical_array(j, P, r, ϑ, φ, k0)
    return [Fr[j]; Fϑ[j]; Fφ[j]]
end

"""
    F_sℓm_cartesian(s,ℓ,m,c,R,k0) -> [Fx;Fy;Fz]

Return normalized vector spherical wave function at R in cartesian coordinates 
"""
function F_sℓm_cartesian(s, ℓ, m, P, R, k0)
    r = norm(R)
    φ = mod(atan(R[2], R[1]), 2pi)
    ϑ = abs(mod(atan(sqrt(R[1]^2 + R[2]^2), R[3]) + pi, 2pi) - pi)
    Fspherical = F_sℓm_spherical(s, ℓ, m, P, r, ϑ, φ, k0)
    sint, cost = sincos(ϑ)
    sinp, cosp = sincos(φ)
    Fcartesian =
        [
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
    Jmaxx::Integer,
    P::PropagationType,
    R::AbstractVector{T},
    k0::T,
) where {T<:Number}
    r = norm(R)
    φ = zero(T)
    ϑ = zero(T)
    if k0 * r > 100 * eps()
        φ = (mod(atan(R[2], R[1]), 2pi))
        ϑ = abs(mod(atan(sqrt(R[1]^2 + R[2]^2), R[3]) + pi, 2pi) - pi)
        # if ϑ > pi
        #     ϑ=(2*pi-ϑ)
        # end 
    end
    Fr, Fϑ, Fφ = F_sℓm_spherical_array(Jmaxx, P, r, ϑ, φ, k0)
    sint, cost = sincos(ϑ)
    sinp, cosp = sincos(φ)
    Fcartesian =
        [Fr Fϑ Fφ] * [
            sint*cosp sint*sinp cost
            cost*cosp cost*sinp -sint
            -sinp cosp 0
        ]
    return convert.(Complex{T}, Fcartesian[:, 1]),
    convert.(Complex{T}, Fcartesian[:, 2]),
    convert.(Complex{T}, Fcartesian[:, 3])
end

function curlF_sℓm_cartesian_array(
    Jmaxx::Integer,
    P::PropagationType,
    R::AbstractVector{T},
    k0::T,
) where {T<:Number}
    Fx, Fy, Fz = F_sℓm_cartesian_array(Jmaxx, P, R, k0)

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
function K_sℓm_thetazero(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    φ::Real,
    outputtype::Type{C} = ComplexF64,
) where {C<:Complex}
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

"""
    K_sℓm(s,ℓ,m,ϑ,φ) -> Kϑ,Kφ

Return farfield vector spherical wave function in spherical coordinates 
"""
function K_sℓm(
    s::Integer,
    ℓ::Integer,
    m::Integer,
    ϑ::Number,
    φ::Real,
    outputtype::Type{C} = ComplexF64,
) where {C<:Complex}
    j = sℓm_to_j(s, ℓ, m)
    Kϑ, Kφ = K_sℓm_array(j, ϑ, φ, outputtype)
    return Kϑ[end], Kφ[end]
end

"""
    K_sℓm_array(Jmaxx,ϑ,φ) -> Kϑ[1:Jmaxx], Kφ[1:Jmaxx]

Return farfield for all vector spherical wave functions up to Jmaxx in spherical coordinates
"""
function K_sℓm_array(
    Jmaxx::Integer,
    ϑ::Number,
    φ::Real,
    outputtype::Type{C} = ComplexF64,
) where {C<:Complex}
    Kϑ, Kφ = K_sℓm_array(Jmaxx::Integer, ϑ::Number, [φ], outputtype)
    return Kϑ[:, 1], Kφ[:, 1]
end
function K_sℓm_array(
    Jmaxx::Integer,
    ϑ::Number,
    φvec::Array{<:Real,1},
    outputtype::Type{C} = ComplexF64,
) where {C<:Complex}
    _, Lmax, __ = j_to_sℓm(Jmaxx)
    Jmax = 2 * Lmax * (Lmax + 2)
    Nφ = length(φvec)

    Kϑ = zeros(outputtype, Jmax, Nφ)
    Kφ = zeros(outputtype, Jmax, Nφ)
    if (abs(ϑ) > 10 * eps()) && (abs(ϑ - pi) > 10 * eps())
        zℓ = [outputtype((1im)^(ℓ + 1) / (sqrt(ℓ * (ℓ + 1)))) for ℓ = 1:Lmax]
        dzℓ = [outputtype((1im)^(ℓ) / (sqrt(ℓ * (ℓ + 1)))) for ℓ = 1:Lmax]

        for m = 1:Lmax
            P1, P2, __ = legendre_deps_array(m, Lmax, ϑ)
            exp_jmφ = outputtype.(cis.(m * φvec))
            exp_minusjmφ = conj(exp_jmφ) * (-1)^m
            for ℓ = m:Lmax
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
        for ℓ = 1:Lmax # m==0
            j = sℓm_to_j(1, ℓ, 0)
            Kϑ[j, :] = 1im * zℓ[ℓ] * P1[ℓ+1] * ones(outputtype, Nφ)
            Kφ[j, :] = -zℓ[ℓ] * P2[ℓ+1] * ones(outputtype, Nφ)

            j = sℓm_to_j(2, ℓ, 0)
            Kϑ[j, :] = dzℓ[ℓ] * P2[ℓ+1] * ones(outputtype, Nφ)
            Kφ[j, :] = 1im * dzℓ[ℓ] * P1[ℓ+1] * ones(outputtype, Nφ)
        end

        return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
    else
        if abs(ϑ) <= 10 * eps()
            for ℓ = 1:Lmax
                j = sℓm_to_j(1, ℓ, -1)
                for (kφ, φ) in enumerate(φvec)
                    Kϑ[j, kφ], Kφ[j, kφ] = K_sℓm_thetazero(1, ℓ, -1, φ)
                    Kϑ[j+1, kφ], Kφ[j+1, kφ] = K_sℓm_thetazero(2, ℓ, -1, φ)

                    Kϑ[j+4, kφ], Kφ[j+4, kφ] = K_sℓm_thetazero(1, ℓ, 1, φ)
                    Kϑ[j+5, kφ], Kφ[j+5, kφ] = K_sℓm_thetazero(2, ℓ, 1, φ)
                end
            end
            return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
        elseif (abs(ϑ - pi) <= 10 * eps())
            for ℓ = 1:Lmax
                j = sℓm_to_j(1, ℓ, -1)
                for (kφ, φ) in enumerate(φvec)
                    Kϑ[j, kφ], Kφ[j, kφ] = K_sℓm_thetapi(1, ℓ, -1, φ)
                    Kϑ[j+1, kφ], Kφ[j+1, kφ] = K_sℓm_thetapi(2, ℓ, -1, φ)

                    Kϑ[j+4, kφ], Kφ[j+4, kφ] = K_sℓm_thetapi(1, ℓ, 1, φ)
                    Kϑ[j+5, kφ], Kφ[j+5, kφ] = K_sℓm_thetapi(2, ℓ, 1, φ)
                end
            end
            return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
        end
        return Kϑ[1:Jmaxx, :], Kφ[1:Jmaxx, :]
    end


end
