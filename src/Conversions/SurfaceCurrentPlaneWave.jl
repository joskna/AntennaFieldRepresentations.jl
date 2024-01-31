function convertrepresentation(
    T::Type{FarfieldPattern{C}}, currents::AbstractSurfaceCurrentDensity, L::Integer, k₀::Number
) where {C<:Complex}
    _, θvec, ϕvec = samplingrule(L)
    Eθ, Eϕ = farfield(currents, θvec, ϕvec, k₀)
    return FarfieldPattern(L, convert.(C, Eθ), convert.(C, Eϕ))
end

function convertrepresentation(T::Type{FarfieldPattern{C}}, currents::AbstractSurfaceCurrentDensity, k₀::Number) where {C<:Complex}
    ϵ = 1e-7

    rmax = 0.0
    for pos in currents.functionspace.geo.vertices
        if norm(pos) > rmax
            rmax = norm(pos)
        end
    end
    rmax *= 2
    L = maximum([3, Int(ceil(k₀ * rmax + 1.8 * (log10(1 / ϵ))^(2 / 3) * (k₀ * rmax)^(1 / 3)))])
    return convertrepresentation(T, currents, L, k₀)
end
