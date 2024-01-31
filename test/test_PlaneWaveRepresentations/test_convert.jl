using LinearAlgebra
c0 = 299792458
f = 1.0e9

λ = c0 / f
k₀ = 2 * pi / λ

Lmax=10
Jmax= sℓm_to_j(2,Lmax,Lmax)
α=RadiatingSphericalExpansion(zeros(ComplexF64,Jmax))

for ℓ in 1:Lmax 
    α.coefficients[sℓm_to_j(1,ℓ,-1)]=complex(0.25)
    α.coefficients[sℓm_to_j(1,ℓ, 1)]=complex(0.25)
    α.coefficients[sℓm_to_j(2,ℓ,-1)]=complex(0.25)
    α.coefficients[sℓm_to_j(2,ℓ, 1)]=-complex(0.25)
end
Lpattern=maximum([2*Lmax, 6])
_,θvec,ϕvec=  samplingrule(Lpattern)

Fθ=zeros(ComplexF64,length(θvec),length(ϕvec))
Fϕ=zeros(ComplexF64,length(θvec),length(ϕvec))

for kθ in eachindex(θvec) 
    for kϕ in  eachindex(ϕvec)
    Fθ[kθ,kϕ], Fϕ[kθ,kϕ]= farfield(α, θvec[kθ], ϕvec[kϕ])
    end
end


storedPattern=FarfieldPattern(Lpattern,Fθ, Fϕ)

convertedPattern= converttype(PlaneWaveSpectrum{elementtype(storedPattern)}, storedPattern)
reconvertedPattern= converttype(FarfieldPattern{elementtype(convertedPattern)}, convertedPattern)

@test all(storedPattern.Eθ==convertedPattern.Eθ)
@test all(storedPattern.Eϕ==convertedPattern.Eϕ)
@test isa(convertedPattern, PlaneWaveSpectrum)


@test all(storedPattern.Eθ==reconvertedPattern.Eθ)
@test all(storedPattern.Eϕ==reconvertedPattern.Eϕ)
@test isa(reconvertedPattern, FarfieldPattern)

revertedPattern=revertdirection(storedPattern)

for k in 1:(storedPattern.L + 1)  # loop over half of ϕ

    # first half of ϕ values
    for kk in 1:(storedPattern.L + 1) # loop over θ
        @test revertedPattern.Eθ[kk, k] == storedPattern.Eθ[end - kk + 1, k + storedPattern.L + 1]
        @test revertedPattern.Eϕ[kk, k] == -storedPattern.Eϕ[end - kk + 1, k + storedPattern.L + 1]
    end

    # second half of ϕ values
    for kk in 1:(storedPattern.L + 1) # loop over θ
        @test revertedPattern.Eθ[kk, k + storedPattern.L + 1] == storedPattern.Eθ[end - kk + 1, k]
        @test revertedPattern.Eϕ[kk, k + storedPattern.L + 1] == -storedPattern.Eϕ[end - kk + 1, k]
    end
end
