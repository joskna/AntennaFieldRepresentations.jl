function convertrepresentation(
    T::Type{RadiatingSphericalExpansion{C}}, dipoles::Array{<:AbstractDipole,1}, k0::Number
) where {C<:Complex}
    ϵ = 1e-7

    rmax = 0.0
    for dipole in dipoles
        if norm(dipole.pos) > rmax
            rmax = norm(dipole.pos)
        end
    end
    rmax *= 2
    L = maximum([3, Int(ceil(k0 * rmax + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k0 * rmax)^(1 / 3)))])
    # Jmax=sℓm_to_j(2,L,L)
    return convertrepresentation(T, dipoles, L, k0)

end
function convertrepresentation(
    T::Type{IncidentSphericalExpansion{C}}, dipoles::Array{<:AbstractDipole,1}, k0::Number
) where {(C <: Complex)}
    ϵ = 1e-7

    rmin = Inf
    for dipole in dipoles
        if norm(dipole.pos) < rmin
            rmin = norm(dipole.pos)
        end
    end
    L = Int(ceil((k0 * rmin + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k0 * rmin)^(1 / 3)) * 2.5))
    # Jmax=sℓm_to_j(2,L,L)
    return convertrepresentation(T, dipoles, L, k0)

end
function convertrepresentation(T::Type{<:AbstractSphericalExpansion}, dipoles::Array{<:AbstractDipole,1}, Lmax::Integer, k0::Number)
    # function convert(T::Type{<:AbstractSphericalExpansion}, dipoles::Array{HertzDipole{F},1}, Jmax::Integer, k0::Real) where{F<:Number}
    Trec = reciprocaltype(T)
    Jmax = sℓm_to_j(2, Lmax, Lmax)

    C = elementtype(T)
    α = zeros(C, Jmax)
    α_hlp = zeros(C, Jmax)
    for dipole in dipoles
        α_hlp .+= _dipole_spherical_innerprod(dipole, Jmax, Trec, k0)
    end
    for (j, coefficient) in enumerate(α_hlp)
        s, ℓ, m = j_to_sℓm(j)
        α[sℓm_to_j(s, ℓ, -m)] = (-1)^(m + 1) * coefficient
    end


    # α=(α_hlp)
    return T(α)
end

global sqrtZF = sqrt(Z₀)
function _dipole_spherical_innerprod(dipole::HertzDipole{<:Number}, Jmax::Integer, T::Type{<:AbstractSphericalExpansion}, k0::Number)
    Fx, Fy, Fz = F_sℓm_cartesian_array(Jmax, T, dipole.pos, k0)
    return [Fx Fy Fz] * dipole.dir * dipole.mag * k0 * sqrtZF
end
function _dipole_spherical_innerprod(
    dipole::FitzgeraldDipole{<:Number}, Jmax::Integer, T::Type{<:AbstractSphericalExpansion}, k0::Number
)
    Fx, Fy, Fz = curlF_sℓm_cartesian_array(Jmax, T, dipole.pos, k0)
    return [Fx Fy Fz] * dipole.dir * dipole.mag * complex(0.0, -k0) / sqrtZF
end


