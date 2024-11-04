using AntennaFieldRepresentations
using LinearAlgebra
using Serialization


include("setup_functions.jl")

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = rotate(
    generate_AUTdips(
        collect(-0.25λ:λ/4:0.25λ),
        collect(-0.25λ:λ/4:0λ),
        collect(-0.5λ:λ/4:0.5λ),
        k0,
    ),
    0.7,
    0.9,
    1.3,
)

filenameswe = joinpath("testdata", "swe.afr")
swe = changerepresentation(SphericalWaveExpansion{Radiated}, dipoles)
serialize(filenameswe, swe)
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




    filenameb = joinpath("testdata", string("b_", Jextraθ, "_", Jextraϕ, ".afr"))
    # b = transmit(swe, fs)
    # serialize(filenameb, b)
    b = deserialize(filenameb)


    filenamebarborder =
        joinpath("testdata", string("barborder_", Jextraθ, "_", Jextraϕ, ".afr"))
    # barborder = transmit(swe, fsarborder)
    # serialize(filenamebarborder, barborder)
    barborder = deserialize(filenamebarborder)





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


    @test norm(ffswe - reshape(b, size(ffswe))) / norm(b) < 3e-14
    @test norm(ffswe - reshape(barborder, size(ffswe))) / norm(ffswe) < 1e-14
    # @test norm(ffdip - ffswe) / norm(ffdip) < 3e-14

    stm = AntennaFieldRepresentations.TransmitMap(swe, fs)
    b2 = stm(sphcoeffs)
    swe .= sphcoeffs

    @test norm(b - b2) / norm(b) < 1e-16

    # if Jextraθ >= 0 && Jextraϕ >= 0
    #     stm_inv = AntennaFieldRepresentations.inverse(stm)
    #     αret = stm_inv(vec(fs.S21values))

    #     @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 3e-14
    # end

    stm_ad = adjoint(stm)

    # filenameA = joinpath("testdata",string("Amat_", Jextraθ, "_", Jextraϕ, ".afr"))
    # # A = zeros(ComplexF64, size(stm))
    # # Threads.@threads for k = 1:length(sphcoeffs)
    # #     stm_tmp=deepcopy(stm)
    # #     x = zeros(ComplexF64, length(sphcoeffs))
    # #     x[k] = 1
    # #     A[:, k] .= stm_tmp * x
    # # end
    # # serialize(filenameA, A)
    # A = deserialize(filenameA)

    # Aᴴ = similar(A')
    # Threads.@threads for k = 1:length(b)
    #     local y = zeros(ComplexF64, length(b))
    #     stm_ad_tmp=deepcopy(stm_ad)
    #     y[k] = 1
    #     Aᴴ[:, k] .= stm_ad_tmp * y
    # end

    # @test norm(A' .- Aᴴ) / norm(A) < 1e-15

    #    stmarborder = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsarborder)
    #    stmarborder_ad =adjoint(stmarborder)

    #    A=zeros(ComplexF64, size(stmarborder))
    #    Aᴴ=similar(A')

    #    for k in 1: length(sphcoeffs)
    #        x= zeros(ComplexF64, length(sphcoeffs))
    #        x[k]=1
    #        A[:,k]= stmarborder * x 
    #    end

    #    for k in 1: length(barborder)
    #        y= zeros(ComplexF64, length(barborder))
    #        y[k]=1
    #        Aᴴ[:,k]= stmarborder_ad * y 
    #    end

    #    @test norm(A'.-Aᴴ)/ norm(A) < 1e-15



    gausssamplingstrategy =
        GaussLegendreθRegularϕSampling(Lmax + 1 + Jextraθ, 2Lmax + 2 + Jextraϕ)
    fsgauss = SphericalFieldSampling(gausssamplingstrategy, αincext)
    # fsarbordergauss = SphericalFieldSampling(gausssamplingstrategy, αinc_arborder) 
    stm.swe .= sphcoeffs
    bgauss = transmit(swe, fsgauss)
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

    # @test norm(ffdip - bgauss) / norm(ffdip) < 3e-14
    @test norm(ffswe - reshape(bgauss, size(ffswe))) / norm(ffswe) < 3e-14
    # @test norm(ffswe - bgaussarborder) / norm(ffswe) < 1e-14

    stmgauss = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsgauss)
    b2 = stmgauss(sphcoeffs)
    swe .= sphcoeffs

    @test norm(bgauss - b2) < 1e-16

    # if Jextraθ >= 0 && Jextraϕ >= 0
    #     # a = AntennaFieldRepresentations.fastsphericalinverse(
    #     #     ffswe,
    #     #     fsgauss.incidentcoefficients,
    #     #     θs,
    #     #     θweights,
    #     # )
    #     # @test norm(a[1:length(sphcoeffs)] - sphcoeffs) / norm(sphcoeffs) < 4e-14

    #     stm_inv = AntennaFieldRepresentations.inverse(stmgauss)
    #     αret = stm_inv(vec(ffswe))
    #     @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 4e-14

    # end

    #     stmgauss_ad =adjoint(stmgauss)
    #     A=zeros(ComplexF64, size(stmgauss))
    #     Aᴴ=similar(A')

    #    for k in 1: length(sphcoeffs)
    #        x= zeros(ComplexF64, length(sphcoeffs))
    #        x[k]=1
    #        A[:,k]= stmgauss * x 
    #    end

    #    for k in 1: length(bgauss)
    #        y= zeros(ComplexF64, length(bgauss))
    #        y[k]=1
    #        Aᴴ[:,k]= stmgauss_ad * y 
    #    end

    #    @test norm(A'.-Aᴴ)/ norm(A) < 1e-15
end
