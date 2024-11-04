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

filenameswe = joinpath("testdata", "swe.afr")
# swe = changerepresentation(SphericalWaveExpansion{Radiated}, dipoles)
# serialize(filenameswe, swe)
swe = deserialize(filenameswe)

sphcoeffs = deepcopy(swe.coefficients)

_, Lmax, __ = j_to_sℓm(length(swe))

αinc = AntennaFieldRepresentations.αinc_planewave(Lmax)
αinc_arborder = SphericalCoefficients(αinc)

swe .= sphcoeffs

αincext = AntennaFieldRepresentations.αinc_planewave(Lmax + 5)

for Jextraθ = -2:3, Jextraϕ = -2:3



    Jθ = 2 * Lmax + 1 + Jextraθ
    Jϕ = 2 * Lmax + 1 + Jextraϕ

    # Jθ= 2*Lmax + 5
    # Jϕ= 2*Lmax + 3

    # Jθ= 10
    # Jϕ= 10
    regularsamplingstrategy = RegularθRegularϕSampling(Jθ, Jϕ)
    fs = SphericalFieldSampling(regularsamplingstrategy, αincext)
    fsarborder = SphericalFieldSampling(regularsamplingstrategy, αinc_arborder)




    # filenameb=joinpath("testdata",string("b_", Jextraθ, "_", Jextraϕ, ".afr"))
    # b = transmit(swe, fs)
    # serialize(filenameb, b)
    # b=deserialize(filenameb)


    # filenamebarborder=joinpath("testdata",string("barborder_", Jextraθ, "_", Jextraϕ, ".afr"))
    # barborder = transmit(swe, fsarborder)
    # serialize(filenamebarborder, barborder)
    # barborder = deserialize(filenamebarborder)





    # @btime b = reshape(transmit(swe, fs), 41,81,2);


    θweights, ϕweights, θs, ϕs =
        AntennaFieldRepresentations.weightsandsamples(regularsamplingstrategy)
    filenameffswe = joinpath("testdata", string("ffswe_", Jextraθ, "_", Jextraϕ, ".afr"))
    # ffswe = zeros(ComplexF64, size(fs.S21values))
    # # ffdip= zeros(ComplexF64, size(fs.S21values))
    # for k in eachindex(θs), kk in eachindex(ϕs)
    #     local Eθ, Eϕ = farfield(swe, (θs[k], ϕs[kk]))
    #     ffswe[k, kk, 1] = Eθ
    #     ffswe[k, kk, 2] = Eϕ

    #     # local Eθ, Eϕ = farfield(dipoles, (θs[k], ϕs[kk]))
    #     # ffdip[k,kk, 1] = Eθ
    #     # ffdip[k,kk, 2] = Eϕ
    # end
    # serialize(filenameffswe, ffswe)
    ffswe = deserialize(filenameffswe)


    # @test norm(ffdip - ffswe) / norm(ffdip) < 3e-14

    stm = AntennaFieldRepresentations.SphericalTransmitMap(swe, fs)
    # b2 = stm(sphcoeffs)
    swe .= sphcoeffs



    if Jextraθ >= 0 && Jextraϕ >= 0
        stm_inv = AntennaFieldRepresentations.inverse(stm)
        αret = stm_inv(vec(ffswe))

        @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 3e-14
    end



    gausssamplingstrategy =
        GaussLegendreθRegularϕSampling(Lmax + 1 + Jextraθ, 2Lmax + 2 + Jextraϕ)
    fsgauss = SphericalFieldSampling(gausssamplingstrategy, αincext)
    # fsarbordergauss = SphericalFieldSampling(gausssamplingstrategy, αinc_arborder) 
    stm.swe .= sphcoeffs
    # bgauss = transmit(swe, fsgauss)
    # bgaussarborder = reshape(transmit(swe, fsarbordergauss), size(fsarbordergauss.S21values))
    θweights, ϕweights, θs, ϕs =
        AntennaFieldRepresentations.weightsandsamples(gausssamplingstrategy)

    filenameffswegauss =
        joinpath("testdata", string("ffswegauss", Jextraθ, "_", Jextraϕ, ".afr"))
    # ffswe = zeros(ComplexF64, size(fsgauss.S21values))
    # # ffdip = zeros(ComplexF64, size(fsgauss.S21values))
    # for k in eachindex(θs), kk in eachindex(ϕs)
    #     local Eθ, Eϕ = farfield(swe, (θs[k], ϕs[kk]))
    #     ffswe[k, kk, 1] = Eθ
    #     ffswe[k, kk, 2] = Eϕ

    #     # local Eθ, Eϕ = farfield(dipoles, (θs[k], ϕs[kk]))
    #     # ffdip[k, kk, 1] = Eθ
    #     # ffdip[k, kk, 2] = Eϕ
    # end
    # serialize(filenameffswegauss, ffswe)
    ffswe = deserialize(filenameffswegauss)

    stmgauss = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsgauss)

    if Jextraθ >= 0 && Jextraϕ >= 0

        stm_inv = AntennaFieldRepresentations.inverse(stmgauss)
        αret = stm_inv(vec(ffswe))
        @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 4e-14

    end
end
