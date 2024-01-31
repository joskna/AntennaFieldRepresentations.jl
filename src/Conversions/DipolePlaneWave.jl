function convertrepresentation(T::Type{FarfieldPattern{C}}, dipoles::Array{<:AbstractDipole,1}, k₀::Number) where {C<:Complex}
    ϵ = 1e-7

    rmax = 0.0
    for dipole in dipoles
        if norm(dipole.pos) > rmax
            rmax = norm(dipole.pos)
        end
    end
    rmax *= 2
    L = maximum([3, Int(ceil(k₀ * rmax + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k₀ * rmax)^(1 / 3)))])
    return convertrepresentation(T, dipoles, L, k₀)
end
function convertrepresentation(
    T::Type{FarfieldPattern{C}}, dipoles::Array{<:AbstractDipole,1}, L::Integer, k₀::Number
) where {C<:Complex}
    _, θvec, ϕvec = samplingrule(L)
    Eθ, Eϕ = farfield(dipoles, θvec, ϕvec, k₀)
    return T(L, convert.(C, Eθ), convert.(C, Eϕ))
end

function convertrepresentation(T::Type{PlaneWaveSpectrum{C}}, dipoles::Array{<:AbstractDipole,1}, k₀::Number) where {C<:Complex}
    αinc = convertrepresentation(IncidentSphericalExpansion{C}, dipoles, k₀)
    pws = convertrepresentation(T, αinc)
    return pws
end

function convertrepresentation(
    T::Type{PlaneWaveSpectrum{C}}, dipoles::Array{<:AbstractDipole,1}, L::Integer, k₀::Number
) where {C<:Complex}
    αinc = convertrepresentation(IncidentSphericalExpansion{C}, dipoles, L, k₀)
    pws = convertrepresentation(T, αinc)
    return pws
end

function _convertrepresentation_and_resample(
    T::Type{FarfieldPattern{C}}, dipole::AbstractDipole, L::Integer, k₀::Number
) where {C<:Complex}
    _, θvec, ϕvec = samplingrule(L)
    Eθ::Matrix{C}, Eϕ::Matrix{C} = farfield([dipole], θvec, ϕvec, k₀)
    return T(L, Eθ, Eϕ)
end
