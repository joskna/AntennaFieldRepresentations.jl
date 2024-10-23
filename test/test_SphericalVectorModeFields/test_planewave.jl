#compare Hansen p.341
# Z₀=376.730313669
Z₀ = AntennaFieldRepresentations.Z₀
c0 = AntennaFieldRepresentations.c₀
f = 1.0e9

λ = c0 / f
k0 = 2 * pi / λ

xvec = -2.0*λ:λ/2:2.0*λ
yvec = -2.0*λ:λ/2:2.0*λ
zvec = -2.0*λ:λ/2:2.0*λ
E0 = 1
Lmax = 40
Jmax = sℓm_to_j(2, Lmax, Lmax)
coeffs = zeros(ComplexF64, Int(Jmax))
for ℓ = 1:Lmax
    Q = sqrt(4 * pi) / (sqrt(Z₀) * k0) * 0.5 * sqrt(2ℓ + 1) * E0 * (1im)^(ℓ + 1)
    coeffs[sℓm_to_j(1, ℓ, -1)] = Q
    coeffs[sℓm_to_j(2, ℓ, -1)] = -Q
    coeffs[sℓm_to_j(1, ℓ, 1)] = Q
    coeffs[sℓm_to_j(2, ℓ, 1)] = Q
end
α = SphericalWaveExpansion(Incident(), coeffs, k0)


for x in xvec, y in yvec, z in zvec
    Einc = [
        E0 * cis(k0 * z)
        complex(0.0)
        complex(0.0)
    ]

    Hinc = [
        complex(0.0)
        -E0 / Z₀ * cis(k0 * z)
        complex(0.0)
    ]


    E = efield(α, [x; y; z])
    H = hfield(α, [x; y; z])
    # println(Einc./E)
    # println(Hinc./H)
    @test E ≈ Einc
    @test H ≈ Hinc
end

# αrad=converttype(RadiatingSphericalExpansion{ComplexF64},α)
# αabs=converttype(AbsorbedSphericalExpansion{ComplexF64},α)
# αinc=converttype(IncidentSphericalExpansion{ComplexF64},α)

# @test αrad.coefficients ==  α.coefficients
# @test αabs.coefficients ==  α.coefficients
# @test αinc.coefficients ==  α.coefficients

# @test typeof(αrad)==RadiatingSphericalExpansion{ComplexF64}
# @test typeof(αabs)==AbsorbedSphericalExpansion{ComplexF64}
# @test typeof(αinc)==IncidentSphericalExpansion{ComplexF64}
