function sphericalHankel1(ℓ, z)
    # Wolfram alpha spherical hankel function
    return sqrt(pi / (2 * z)) * hankelh1(ℓ + 1//2, z)
end

function oneoverz_deriv_sphericalHankel1(ℓ, z)
    # Hansen Apendix page 316
    return ((ℓ + 1) * sphericalHankel1(ℓ - 1, z) - ℓ * sphericalHankel1(ℓ + 1, z)) / (2 * ℓ + 1)
end

function sphericalHankel2(ℓ, z)
    # Wolfram alpha spherical hankel function
    return sqrt(pi / (2 * z)) * hankelh2(ℓ + 1//2, z)
end

function oneoverz_deriv_sphericalHankel2(ℓ, z)
    # Hansen Apendix page 316
    return ((ℓ + 1) * sphericalHankel2(ℓ - 1, z) - ℓ * sphericalHankel2(ℓ + 1, z)) / (2 * ℓ + 1)
end

function sphericalBessel1(ℓ, z)
    # Wolfram alpha spherical hankel function
    return sqrt(pi / (2 * z)) * besselj(ℓ + 1//2, z)
end

function oneoverz_deriv_sphericalBessel1(ℓ, z)
    # Hansen Apendix page 316
    return ((ℓ + 1) * sphericalBessel1(ℓ - 1, z) - ℓ * sphericalBessel1(ℓ + 1, z)) / (2 * ℓ + 1)
end

function zc_ℓ(T::Type{<:AbstractSphericalExpansion}, ℓ, kA)
    return complex(Inf)
end
function zc_ℓ(T::Type{RadiatingSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, sphericalHankel2(ℓ, kA))
end
function zc_ℓ(T::Type{AbsorbedSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, sphericalHankel1(ℓ, kA))
end
function zc_ℓ(T::Type{IncidentSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, sphericalBessel1(ℓ, kA))
end

function oneoverkA_deriv_zc_ℓ(T::Type{<:AbstractSphericalExpansion}, ℓ, kA)
    return complex(Inf)
end
function oneoverkA_deriv_zc_ℓ(T::Type{IncidentSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, oneoverz_deriv_sphericalBessel1(ℓ, kA))
end
function oneoverkA_deriv_zc_ℓ(T::Type{RadiatingSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, oneoverz_deriv_sphericalHankel2(ℓ, kA))
end
function oneoverkA_deriv_zc_ℓ(T::Type{AbsorbedSphericalExpansion{C}}, ℓ, kA) where {C<:Complex}
    return convert(C, oneoverz_deriv_sphericalHankel1(ℓ, kA))
end

function start_radial_recurrence(T::Type{RadiatingSphericalExpansion{C}}, kA) where {C<:Complex}
    expfac = cis(-kA)
    zℓ1 = 1im * expfac / kA
    zℓ2 = (1im - kA) * expfac / (kA^2)
    return convert(C, zℓ1), convert(C, zℓ2)
end
function start_radial_recurrence(T::Type{IncidentSphericalExpansion{C}}, kA) where {C<:Complex}
    zℓ1 = sin(kA) / kA
    zℓ2 = sin(kA) / (kA^2) - cos(kA) / kA
    return convert(C, zℓ1), convert(C, zℓ2)
end
function start_radial_recurrence(T::Type{AbsorbedSphericalExpansion{C}}, kA) where {C<:Complex}
    expfac = cis(kA)
    zℓ1 = -1im * expfac / kA
    zℓ2 = -expfac * (kA + 1im) / (kA^2)
    return convert(C, zℓ1), convert(C, zℓ2)
end
function start_radial_recurrence(T::Type{UnorthodoxSphericalExpansion{C}}, kA) where {C<:Complex}
    zℓ1 = -cos(kA) / kA
    zℓ2 = -cos(kA) / (kA^2) - sin(kA) / kA
    return convert(C, zℓ1), convert(C, zℓ2)
end

function R_dependencies_array(T::Type{<:AbstractSphericalExpansion}, Lmax, kA)
    zℓ = zeros(elementtype(T), Lmax + 2)
    dzℓ = zeros(elementtype(T), Lmax + 1)

    zℓ[1], zℓ[2] = start_radial_recurrence(T, kA)
    dzℓ[1] = oneoverkA_deriv_zc_ℓ(T, 0, kA)

    if Lmax > 1
        for ℓ in 2:Lmax
            zℓ[ℓ + 1] = (2 * ℓ - 1) / kA * zℓ[ℓ] - zℓ[ℓ - 1]
        end

        for ℓ in 2:Lmax
            dzℓ[ℓ] = (((ℓ) * zℓ[ℓ - 1] - (ℓ - 1) * zℓ[ℓ + 1]) / (2 * ℓ - 1))
        end
    else
        for ℓ in 1:Lmax
            zℓ[ℓ] = zc_ℓ(T, ℓ - 1, kA)
            dzℓ[ℓ] = oneoverkA_deriv_zc_ℓ(T, ℓ - 1, kA)
        end
    end
    return zℓ[1:(end - 1)], dzℓ
end