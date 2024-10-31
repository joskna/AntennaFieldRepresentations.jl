using AntennaFieldRepresentations
using LinearAlgebra
using Serialization


# include("setup_functions.jl")

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

# dipoles = rotate(
#     generate_AUTdips(
#         collect(-0.25λ:λ/4:0.25λ),
#         collect(-0.25λ:λ/4:0λ),
#         collect(-0.5λ:λ/4:0.5λ),
#         k0,
#     ),
#     0.7,
#     0.9,
#     1.3,
# )

filenameswe=joinpath("testdata","swe.afr")
# swe = changerepresentation(SphericalWaveExpansion{Radiated}, dipoles)
# serialize(filenameswe, swe)
swe= deserialize(filenameswe)

sphcoeffs = deepcopy(swe.coefficients)

_, Lmax, __ = j_to_sℓm(length(swe))

αinc = AntennaFieldRepresentations.αinc_planewave(Lmax)
αinc_arborder = SphericalCoefficients(αinc)

swe .= sphcoeffs

αincext = AntennaFieldRepresentations.αinc_planewave(Lmax + 5)

for Jextraθ = 0:1, Jextraϕ = 0:1



    Jθ = 2 * Lmax + 1 + Jextraθ
    Jϕ = 2 * Lmax + 1 + Jextraϕ

    # Jθ= 2*Lmax + 5
    # Jϕ= 2*Lmax + 3

    # Jθ= 10
    # Jϕ= 10
    regularsamplingstrategy = RegularθRegularϕSampling(Jθ, Jϕ)
    fs = SphericalFieldSampling(regularsamplingstrategy, αincext)

    stm = AntennaFieldRepresentations.SphericalTransmitMap(swe, fs)
    stm_ad = adjoint(stm)

    filenameA = joinpath("testdata",string("Amat_", Jextraθ, "_", Jextraϕ, ".afr"))
    # A = zeros(ComplexF64, size(stm))
    # Threads.@threads for k = 1:length(sphcoeffs)
    #     stm_tmp=deepcopy(stm)
    #     x = zeros(ComplexF64, length(sphcoeffs))
    #     x[k] = 1
    #     A[:, k] .= stm_tmp * x
    # end
    # serialize(filenameA, A)
    A = deserialize(filenameA)

    a,b = size(A)

    Aᴴ = similar(A')
    Threads.@threads for k = 1:a
        local y = zeros(ComplexF64, a)
        stm_ad_tmp=deepcopy(stm_ad)
        y[k] = 1
        Aᴴ[:, k] .= stm_ad_tmp * y
    end

    @test norm(A' .- Aᴴ) / norm(A) < 1e-15

    fsarborder = SphericalFieldSampling(regularsamplingstrategy, αinc_arborder)
    stmarborder = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsarborder)
    stmarborder_ad = adjoint(stm)
    filenameA = joinpath("testdata",string("Amatarborder_", Jextraθ, "_", Jextraϕ, ".afr"))
    # A = zeros(ComplexF64, size(stmarborder))
    # Threads.@threads for k = 1:length(sphcoeffs)
    #     stm_tmp=deepcopy(stmarborder)
    #     x = zeros(ComplexF64, length(sphcoeffs))
    #     x[k] = 1
    #     A[:, k] .= stm_tmp * x
    # end
    # serialize(filenameA, A)
    A = deserialize(filenameA)

    a,b = size(A)

    Aᴴ = similar(A')
    Threads.@threads for k = 1:a
        local y = zeros(ComplexF64, a)
        stm_ad_tmp=deepcopy(stmarborder_ad)
        y[k] = 1
        Aᴴ[:, k] .= stm_ad_tmp * y
    end

    @test norm(A' .- Aᴴ) / norm(A) < 1e-15

  


    gausssamplingstrategy =
        GaussLegendreθRegularϕSampling(Lmax + 1 + Jextraθ, 2Lmax + 2 + Jextraϕ)
    fsgauss = SphericalFieldSampling(gausssamplingstrategy, αincext)


    stmgauss = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsgauss)
    stmgauss_ad =adjoint(stmgauss)
    filenameA = joinpath("testdata",string("Amatgauss", Jextraθ, "_", Jextraϕ, ".afr"))
    # A = zeros(ComplexF64, size(stmgauss))
    # Threads.@threads for k = 1:length(sphcoeffs)
    #     stm_tmp=deepcopy(stmgauss)
    #     x = zeros(ComplexF64, length(sphcoeffs))
    #     x[k] = 1
    #     A[:, k] .= stm_tmp * x
    # end
    # serialize(filenameA, A)
    A = deserialize(filenameA)

    a,b = size(A)


    Aᴴ = similar(A')
    Threads.@threads for k = 1:a
        local y = zeros(ComplexF64, a)
        stm_ad_tmp=deepcopy(stmgauss_ad)
        y[k] = 1
        Aᴴ[:, k] .= stm_ad_tmp * y
    end

    @test norm(A' .- Aᴴ) / norm(A) < 1e-15

end
