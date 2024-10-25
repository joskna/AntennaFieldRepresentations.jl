using AntennaFieldRepresentations
using LinearAlgebra


function generate_AUTdips(
    xvec::Array{Float64,1},
    yvec::Array{Float64,1},
    zvec::Array{Float64,1},
    k0::Float64,
)

    nx = length(xvec)
    ny = length(yvec)
    nz = length(zvec)

    ycenter = (maximum(yvec) + minimum(yvec)) / 2
    ysize = maximum([abs(maximum(yvec) - ycenter), abs(minimum(yvec) - ycenter)])

    ndips = nx * ny * nz
    positions = Vector{Vector{Float64}}(undef, ndips)
    magnitudes = Vector{ComplexF64}(undef, ndips)

    # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
                # println(index)
                positions[index] = [xvec[kkk], yvec[kk], zvec[k]]
                magnitudes[index] = mag
                # println(pos)


            end
        end
    end
    return HertzArray(
        positions,
        [complex.([0.0, 0.0, 1.0]) for k = 1:length(positions)],
        magnitudes,
        k0,
    )

end

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = rotate(
    generate_AUTdips(
        collect(-0.5λ:λ/4:0.5λ),
        collect(-0.5λ:λ/4:0λ),
        collect(-1λ:λ/4:1λ),
        k0,
    ),
    0.7,
    0.9,
    1.3,
)


swe = changerepresentation(SphericalWaveExpansion{Radiated}, dipoles)
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

    b = transmit(swe, fs)
    barborder = transmit(swe, fsarborder)


    # @btime b = reshape(transmit(swe, fs), 41,81,2);


    θweights, ϕweights, θs, ϕs =
        AntennaFieldRepresentations.weightsandsamples(regularsamplingstrategy)

    ffswe = zeros(ComplexF64, size(fs.S21values))
    # ffdip= zeros(ComplexF64, size(fs.S21values))
    for k in eachindex(θs), kk in eachindex(ϕs)
        local Eθ, Eϕ = farfield(swe, (θs[k], ϕs[kk]))
        ffswe[k, kk, 1] = Eθ
        ffswe[k, kk, 2] = Eϕ

        # local Eθ, Eϕ = farfield(dipoles, (θs[k], ϕs[kk]))
        # ffdip[k,kk, 1] = Eθ
        # ffdip[k,kk, 2] = Eϕ
    end
    @test norm(ffswe - reshape(b, size(ffswe))) / norm(b) < 3e-14
    @test norm(ffswe - reshape(barborder, size(ffswe))) / norm(ffswe) < 1e-14
    # @test norm(ffdip - ffswe) / norm(ffdip) < 3e-14

    stm = AntennaFieldRepresentations.SphericalTransmitMap(swe, fs)
    b2 = stm(sphcoeffs)
    swe .= sphcoeffs

    @test norm(b - b2) / norm(b) < 1e-16

    if Jextraθ >= 0 && Jextraϕ >= 0
        stm_inv = AntennaFieldRepresentations.inverse(stm)
        αret = stm_inv(vec(fs.S21values))

        @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 3e-14
    end



    gausssamplingstrategy =
        GaussLegendreθRegularϕSampling(Lmax + 1 + Jextraθ, 2Lmax + 2 + Jextraϕ)
    fsgauss = SphericalFieldSampling(gausssamplingstrategy, αincext)
    # fsarbordergauss = SphericalFieldSampling(gausssamplingstrategy, αinc_arborder) 
    bgauss = transmit(swe, fsgauss)
    # bgaussarborder = reshape(transmit(swe, fsarbordergauss), size(fsarbordergauss.S21values))
    θweights, ϕweights, θs, ϕs =
        AntennaFieldRepresentations.weightsandsamples(gausssamplingstrategy)

    ffswe = zeros(ComplexF64, size(fsgauss.S21values))
    # ffdip = zeros(ComplexF64, size(fsgauss.S21values))
    for k in eachindex(θs), kk in eachindex(ϕs)
        local Eθ, Eϕ = farfield(swe, (θs[k], ϕs[kk]))
        ffswe[k, kk, 1] = Eθ
        ffswe[k, kk, 2] = Eϕ

        # local Eθ, Eϕ = farfield(dipoles, (θs[k], ϕs[kk]))
        # ffdip[k, kk, 1] = Eθ
        # ffdip[k, kk, 2] = Eϕ
    end
    # @test norm(ffdip - bgauss) / norm(ffdip) < 3e-14
    @test norm(ffswe - reshape(bgauss, size(ffswe))) / norm(ffswe) < 1e-14
    # @test norm(ffswe - bgaussarborder) / norm(ffswe) < 1e-14

    stmgauss = AntennaFieldRepresentations.SphericalTransmitMap(swe, fsgauss)
    b2 = stmgauss(sphcoeffs)
    swe .= sphcoeffs

    @test norm(bgauss - b2) < 1e-16

    if Jextraθ >= 0 && Jextraϕ >= 0
        a = AntennaFieldRepresentations.fastsphericalinverse(
            ffswe,
            fsgauss.incidentcoefficients,
            θs,
            θweights,
        )
        @test norm(a[1:length(sphcoeffs)] - sphcoeffs) / norm(sphcoeffs) < 4e-14

        stm_inv = AntennaFieldRepresentations.inverse(stmgauss)
        αret = stm_inv(vec(ffswe))
        @test norm(αret - sphcoeffs) / norm(sphcoeffs) < 4e-14
    end
end
