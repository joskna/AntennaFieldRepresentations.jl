function rotate(α::AbstractSphericalExpansion, χ::Real, θ::Real, ϕ::Real)
    # coefficients = deepcopy(α.coefficients)
    α_rotated = deepcopy(α)
    _, Lmax, __ = j_to_sℓm(length(α.coefficients))
    Jmax = 2 * Lmax * (Lmax + 2)

    α_rotated.coefficients = similar(α.coefficients, Jmax)
    fill!(α_rotated.coefficients, 0.0)

    for ℓ in 1:Lmax
        d = wignerd(ℓ, -θ)
        for m in (-ℓ):ℓ
            expϕ = cis(-m * ϕ)
            for s in 1:2
                jj = sℓm_to_j(s, ℓ, m)
                for μ in (-ℓ):ℓ
                    j = sℓm_to_j(s, ℓ, μ)
                    if j <= length(α.coefficients)

                        α_rotated.coefficients[jj] += expϕ * d[μ + ℓ + 1, m + ℓ + 1] * cis(-μ * χ) * α.coefficients[j]
                    end
                end
            end
        end
    end

    return α_rotated
end




function αtoβ(α::Array{<:Complex,1})
    β = similar(α, length(α))
    for k in 1:length(α)
        s, ℓ, m = j_to_sℓm(k)
        β[k] = (-1)^(m) * (α[sℓm_to_j(s, ℓ, -m)])
    end
    return β * 0.5
end

function βtoα(β::Array{<:Complex,1})
    α = similar(β, length(β))
    for k in 1:length(β)
        s, ℓ, m = j_to_sℓm(k)
        α[k] = (-1)^(m) * (β[sℓm_to_j(s, ℓ, -m)])
    end
    return α * 2
end
