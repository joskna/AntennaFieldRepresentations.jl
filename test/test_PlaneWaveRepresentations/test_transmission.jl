using Random

c0 = 299792458
f = 1.0e9

λ = c0 / f
k0 = 2 * pi / λ

rng = MersenneTwister(1234);

L=15

for k in 1:3
    αradiated=RadiatingSphericalExpansion(randn(ComplexF64, sℓm_to_j(2,L,L)))
    αincident=IncidentSphericalExpansion(randn(ComplexF64, sℓm_to_j(2,L,L)))

    ffrad=convertrepresentation(FarfieldPattern{ComplexF64}, αradiated, k0)
    pwsinc=convertrepresentation(PlaneWaveSpectrum{ComplexF64}, αincident, k0)

    b1=transmission(αradiated, αincident,k0)
    b2=transmission(pwsinc, ffrad)

    @test b1≈b2
end