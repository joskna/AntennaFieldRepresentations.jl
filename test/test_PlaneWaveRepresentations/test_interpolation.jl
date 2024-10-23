using LinearAlgebra
c0 = 299792458
f = 1.0e9

λ = c0 / f
k₀ = 2 * pi / λ

Lmax = 10
Jmax = sℓm_to_j(2, Lmax, Lmax)
α = RadiatingSphericalExpansion(zeros(ComplexF64, Jmax))

for ℓ = 1:Lmax
    α.coefficients[sℓm_to_j(1, ℓ, -1)] = complex(0.25)
    α.coefficients[sℓm_to_j(1, ℓ, 1)] = complex(0.25)
    α.coefficients[sℓm_to_j(2, ℓ, -1)] = complex(0.25)
    α.coefficients[sℓm_to_j(2, ℓ, 1)] = -complex(0.25)
end
Lpattern = maximum([2 * Lmax, 6])
_, θvec, ϕvec = samplingrule(Lpattern)

Fθ = zeros(ComplexF64, length(θvec), length(ϕvec))
Fϕ = zeros(ComplexF64, length(θvec), length(ϕvec))

for kθ in eachindex(θvec)
    for kϕ in eachindex(ϕvec)
        Fθ[kθ, kϕ], Fϕ[kθ, kϕ] = farfield(α, θvec[kθ], ϕvec[kϕ])
    end
end


storedPattern = FarfieldPattern(Lpattern, Fθ, Fϕ)
orderθ, orderϕ = 14, 14

Lpatternnew = Lpattern + 20
interpattern = resample(
    storedPattern,
    AntennaFieldRepresentations.LocalInterpolation{
        Lpattern,
        Lpatternnew,
        orderθ,
        orderϕ,
        Float64,
    }(),
)
_, θvec, ϕvec = samplingrule(Lpatternnew)

Fθ = zeros(ComplexF64, length(θvec), length(ϕvec))
Fϕ = zeros(ComplexF64, length(θvec), length(ϕvec))

for kθ in eachindex(θvec)
    for kϕ in eachindex(ϕvec)
        Fθ[kθ, kϕ], Fϕ[kθ, kϕ] = farfield(α, θvec[kθ], ϕvec[kϕ])
    end
end
largepattern = FarfieldPattern(Lpattern, Fθ, Fϕ)

@test interpattern.Eθ ≈ largepattern.Eθ rtol = 5e-4
@test interpattern.Eϕ ≈ largepattern.Eϕ rtol = 5e-4


