"""
αinc_PW(L::Integer)


Return incident spherical mode coefficients up to mode order L
of an x-polarized plane wave traveling in negative z-direction.
Compare Hansen: "Spherical Near-Field Measurements" Appendix A1.6 
"""
function αinc_PW(L::Integer)
    αin = zeros(ComplexF64, sℓm_to_j(2, L, L))
    for ℓ in 1:L

        αin[sℓm_to_j(1, ℓ, 1)] = -(1im)^(ℓ) * sqrt((2 * ℓ + 1) / pi)
        αin[sℓm_to_j(2, ℓ, 1)] = -(1im)^(ℓ) * sqrt((2 * ℓ + 1) / pi)
        αin[sℓm_to_j(1, ℓ, -1)] = -(1im)^(ℓ) * sqrt((2 * ℓ + 1) / pi)
        αin[sℓm_to_j(2, ℓ, -1)] = (1im)^(ℓ) * sqrt((2 * ℓ + 1) / pi)
    end
    return αin * sqrt(Z₀) * 0.5
end


function _convertrepresentation_and_resample(
    T::Type{FarfieldPattern{C2}}, α::RadiatingSphericalExpansion{C1}, L::Integer
) where {C1<:Complex,C2<:Complex}
    _, Lmax, __ = j_to_sℓm(length(α.coefficients))
    β_aut = zeros(C1, sℓm_to_j(2, Lmax, Lmax)) # ensure that β_aut has correct length
    β_aut[1:length(α.coefficients)] = αtoβ(α.coefficients)
    _, θvec, __ = samplingrule(L)
    S12 = Wackerforward(αinc_PW(Lmax), β_aut, θvec, 2 * L + 2)
    return T(L, convert.(C2, S12[:, :, 1]), convert.(C2, -S12[:, :, 2]))
end
function _convertrepresentation_and_resample(
    T::Type{FarfieldPattern{C2}}, α::RadiatingSphericalExpansion{C1}, L::Integer, ::Real
) where {C1<:Complex,C2<:Complex}
return _convertrepresentation_and_resample(T, α, L)
end
function _convertrepresentation_and_resample(
    T::Type{PlaneWaveSpectrum{C2}}, α::IncidentSphericalExpansion{C1}, L::Integer
) where {C1<:Complex,C2<:Complex}
    pws = _convertrepresentation_and_resample(FarfieldPattern{C2}, converttype(RadiatingSphericalExpansion{C1}, α), L)
    return converttype(T, pws)
end
function convertrepresentation(T::Type{FarfieldPattern{C2}}, α::RadiatingSphericalExpansion{C1}) where {C1<:Complex,C2<:Complex}
    _, Lmax, __ = j_to_sℓm(length(α.coefficients))
    β_aut = zeros(C1, sℓm_to_j(2, Lmax, Lmax)) # ensure that β_aut has correct length
    β_aut[1:length(α.coefficients)] = αtoβ(α.coefficients)
    _, θvec, __ = samplingrule(Lmax)
    S12 = Wackerforward(αinc_PW(Lmax), β_aut, θvec, 2 * Lmax + 2)
    return T(Lmax, convert.(C2, S12[:, :, 1]), convert.(C2, -S12[:, :, 2]))
end
function convertrepresentation(
    T::Type{FarfieldPattern{C2}}, α::RadiatingSphericalExpansion{C1}, ::Real
) where {C1<:Complex,C2<:Complex}
    return convertrepresentation(T, α)
end

function convertrepresentation(T::Type{RadiatingSphericalExpansion{C}}, FF::FarfieldPattern) where {C<:Complex}
    a, b = size(FF.Eθ)

    S12 = zeros(ComplexF64, a, b, 2)
    S12[:, :, 1] = FF.Eθ
    S12[:, :, 2] = -FF.Eϕ
    θweights, θvec, __ = samplingrule(FF.L)
    return T(βtoα(Wacker(S12, αinc_PW(FF.L), θvec, θweights)))
end
function convertrepresentation(T::Type{RadiatingSphericalExpansion{C}}, FF::FarfieldPattern, ::Float64) where {C<:Complex}
    return convertrepresentation(T, FF)
end

function convertrepresentation(T::Type{IncidentSphericalExpansion{C}}, PWS::PlaneWaveSpectrum) where {C<:Complex}
    α = convertrepresentation(RadiatingSphericalExpansion{C}, converttype(FarfieldPattern{C}, PWS))
    return converttype(T, α)
end

function convertrepresentation(T::Type{PlaneWaveSpectrum{C2}}, α::IncidentSphericalExpansion{C1}, k0) where {C1<:Complex,C2<:Complex}
    pws = convertrepresentation(FarfieldPattern{C2}, converttype(RadiatingSphericalExpansion{C1}, α), k0)
    return converttype(T, pws)
end

function convertrepresentation(T::Type{IncidentSphericalExpansion{C}}, PW::PlaneWave, Lmax::Integer) where {C<:Complex}
    # calculate theta and phi angle of incident plane wave
    θ = atan(sqrt(PW.kvec[1]^2 + PW.kvec[2]^2), PW.kvec[3])# result is always positive
    ϕ = mod2pi(atan(PW.kvec[2], PW.kvec[1]))
    sint, cost = sincos(θ)
    sinp, cosp = sincos(ϕ)
    # eᵣ = [cosp * sint; sinp * sint; cost]
    eθ = [cosp * cost; sinp * cost; -sint]
    eϕ = [-sinp; cosp; 0]
    Kϑ, Kφ = K_sℓm_array(sℓm_to_j(2, Lmax, Lmax), θ, ϕ)
    # Kϑ=αtoβ(Kϑ)
    # Kφ=αtoβ(Kφ)
    coeffs = (conj.(Kϑ) * udot(eθ, PW.pol) + conj.(Kφ) * udot(eϕ, PW.pol)) * PW.mag * complex(0.0, 4 * pi / (norm(PW.kvec) * sqrt(Z₀)))

    return T(convert.(C, coeffs))
end



