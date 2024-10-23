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

rotatedα = rotate(α, π / 3, π / 5, π / 7)

Fθ = zeros(ComplexF64, length(θvec), length(ϕvec))
Fϕ = zeros(ComplexF64, length(θvec), length(ϕvec))

for kθ in eachindex(θvec)
    for kϕ in eachindex(ϕvec)
        Fθ[kθ, kϕ], Fϕ[kθ, kϕ] = farfield(rotatedα, θvec[kθ], ϕvec[kϕ])
    end
end

rotatedstoredPattern = FarfieldPattern(Lpattern, Fθ, Fϕ)

rotatedPattern = rotate(storedPattern, π / 3, π / 5, π / 7)
@test all(isapprox.(rotatedPattern.Eθ, rotatedstoredPattern.Eθ, rtol = 1e-1))
@test all(isapprox.(rotatedPattern.Eϕ, rotatedstoredPattern.Eϕ, rtol = 1e-1))

@test norm(rotatedstoredPattern.Eθ - rotatedPattern.Eθ) / norm(rotatedstoredPattern.Eθ) <
      1e-3
@test norm(rotatedstoredPattern.Eϕ - rotatedPattern.Eϕ) / norm(rotatedstoredPattern.Eϕ) <
      1e-3

rotatedPattern = rotate(storedPattern, π / 3, π / 5, π / 7, 30, 30)
@test all(isapprox.(rotatedPattern.Eθ, rotatedstoredPattern.Eθ, rtol = 9e-5))
@test all(isapprox.(rotatedPattern.Eϕ, rotatedstoredPattern.Eϕ, rtol = 9e-5))

@test norm(rotatedstoredPattern.Eθ - rotatedPattern.Eθ) / norm(rotatedstoredPattern.Eθ) <
      1e-6
@test norm(rotatedstoredPattern.Eϕ - rotatedPattern.Eϕ) / norm(rotatedstoredPattern.Eϕ) <
      1e-6
