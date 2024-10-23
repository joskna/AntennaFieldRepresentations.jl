using LinearAlgebra


# copy from ElectromagneticDipoleUtils to not depend on this library
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

    dipoles = Array{HertzDipole{Float64},1}(undef, ndips)# Array(HertzDipole{Float64},ndips)
    # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
                # println(index)
                pos = [xvec[kkk], yvec[kk], zvec[k]]
                # println(pos)
                dipoles[index] =
                    HertzDipole(pos, [0.0 + 0.0im, 0.0 + 0.0im, 1.0 + 0.0im], mag / ndips)

            end
        end
    end
    return dipoles

end
# copy from ElectromagneticDipoleUtils to not depend on this library
function dipoles_on_sphere(
    deltatheta::Float64,
    deltaphi::Float64,
    radius::T;
    tangential = true,
    pole = false,
) where {T<:Number}
    phi = 0:deltaphi:2*pi-deltaphi
    theta = deltatheta/2:deltatheta:pi-deltatheta/2
    if pole == true
        theta = 0:deltatheta:pi
    end


    numtheta = length(theta)
    numphi = length(phi)

    fac = 1 / (numtheta * numphi)

    dipoles = Array{HertzDipole{typeof(radius)},1}(undef, 2 * numtheta * numphi)
    # two tangential dipoles at each location

    if tangential == false
        # return also  radial dipoles  
        dipoles = Array{HertzDipole{typeof(radius)},1}(undef, 3 * numtheta * numphi)
    end



    for k = 1:numtheta
        cost = cos(theta[k])
        sint = sin(theta[k])

        for kk = 1:numphi
            cosp = cos(phi[kk])
            sinp = sin(phi[kk])

            index = numphi * (k - 1) + kk
            x = radius * sint * cosp
            y = radius * sint * sinp
            z = radius * cost


            e_θ = [cost * cosp; cost * sinp; -sint]
            e_ϕ = [-sinp; cosp; 0.0]
            dipoles[index] = HertzDipole([x; y; z], complex(e_θ), 1.0 + 0.0im)
            dipoles[index+numtheta*numphi] =
                HertzDipole([x; y; z], complex(e_ϕ), fac + 0.0im)
            if tangential == false
                e_r = [sint * cosp; sint * sinp; cost]
                dipoles[index+2*numtheta*numphi] =
                    HertzDipole([x; y; z], complex(e_r), fac + 0.0im)
            end
        end

    end

    return dipoles

end

# using ElectromagneticUtils
# ## 1GHz
c0 = 299792458
f = 1.0e9

λ = c0 / f
k0 = 2 * pi / λ
# println(k0)

@testset verbose = true "Conversions Dipole -> RadiatingSpherical" begin
    rad_aut = 0.7 * λ

    rad_meas = 1.51 * rad_aut

    x_aut = sqrt(rad_aut^2 - (λ / 4)^2) / sqrt(2)
    y_aut = x_aut * 1.25



    autdips = generate_AUTdips(
        collect(-λ/4:λ/3:λ/4),
        collect(-x_aut/2:λ/3:x_aut/2),
        collect(-y_aut/2:λ/3:y_aut/2),
        k0,
    )

    α = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, autdips, k0)
    θvec = [0; 2 / 3 * pi]
    ϕvec = collect(0:135:359) / 180 * π
    Fθ, Fϕ = farfield(autdips, θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(α, θvec, ϕvec)


    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))

    xtest = -1.25*x_aut:0.75*λ:1.25*x_aut
    ytest = -1.25*y_aut:0.75:1.25*y_aut
    ztest = [rad_meas]
    for x in xtest, y in ytest, z in ztest
        R = [x, y, z]
        E1 = efield(autdips, R, k0)
        H1 = hfield(autdips, R, k0)
        E2, H2 = ehfield(α, R, k0)
        @test all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
        @test all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))
    end

    ztest = -1.25*x_aut:0.75*λ:1.25*x_aut
    ytest = -1.25*y_aut:0.75*λ:1.25*y_aut
    xtest = [-rad_meas; rad_meas]
    for x in xtest, y in ytest, z in ztest
        R = [x, y, z]
        E1 = efield(autdips, R, k0)
        H1 = hfield(autdips, R, k0)
        E2, H2 = ehfield(α, R, k0)
        @test all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
        @test all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))
    end

    # xtest=-1.25*x_aut:λ/2:1.25*x_aut
    # ztest=-1.25*y_aut:λ/2:1.25*y_aut
    # ytest=[-rad_meas;rad_meas]
    # for x in xtest, y in ytest, z in ztest
    #     R=[x,y,z]
    #     E1= efield(autdips, R, k0)
    #     H1= hfield(autdips, R, k0)
    #     E2, H2=ehfield(α, R, k0)
    #     @test all(isapprox.(E1,E2; rtol=2e-4, atol=1e-6))
    #     @test all(isapprox.(H1,H2; rtol=2e-4, atol=1e-6))

    #     if !all(isapprox.(E1,E2; rtol=2e-4, atol=1e-6))
    #         println(R)
    #         println("E1:", E1)
    #         println("E2:", E2)
    #     end
    #     if !all(isapprox.(H1,H2; rtol=1e-4, atol=1e-6))
    #         println(R)
    #         println("H1:", H1)
    #         println("H2:", H2)
    #     end
    # end

    autdipsF = converttype.(FitzgeraldDipole{Float64}, autdips)
    αF = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, autdipsF, k0)
    Fθ, Fϕ = farfield(autdipsF, θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αF, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))


    xtest = -1.25*x_aut:0.75*λ:1.25*x_aut
    ztest = -1.25*y_aut:0.75*λ:1.25*y_aut
    ytest = [-rad_meas; rad_meas]
    for x in xtest, y in ytest, z in ztest
        R = [x, y, z]
        E1 = efield(autdipsF, R, k0)
        H1 = hfield(autdipsF, R, k0)
        E2, H2 = ehfield(αF, R, k0)
        @test all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
        @test all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))

        if !all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
            println(R)
            println("E1:", E1)
            println("E2:", E2)
        end
        if !all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))
            println(R)
            println("H1:", H1)
            println("H2:", H2)
        end
    end

    autdipsHF = [autdips; autdipsF]
    αHF = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, autdipsHF, k0)
    Fθ, Fϕ = farfield(autdipsHF, θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αHF, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))






    dipole = HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(1.0))
    αHertz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αHertz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))


    dipole = HertzDipole([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], complex(1.0))
    αHertz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αHertz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))





    dipole = HertzDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], complex(1.0))
    αHertz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αHertz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))



    dipole = converttype(
        FitzgeraldDipole{Float64},
        HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(1.0)),
    )
    αHertz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αHertz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))


    dipole = converttype(
        FitzgeraldDipole{Float64},
        HertzDipole([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], complex(1.0)),
    )
    αFitz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αFitz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))

    dipole = converttype(
        FitzgeraldDipole{Float64},
        HertzDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], complex(1.0)),
    )
    αFitz = convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, [dipole], k0)

    Fθ, Fϕ = farfield([dipole], θvec, ϕvec, k0)
    Fθ2, Fϕ2 = farfield(αFitz, θvec, ϕvec)

    @test all(isapprox.(Fθ, Fθ2; rtol = 1e-6, atol = 1e-6))
    @test all(isapprox.(Fϕ, Fϕ2; rtol = 1e-6, atol = 1e-6))
end

@testset verbose = true "Conversions Dipole -> IncidentSpherical" begin
    Hertzobsdips = dipoles_on_sphere(45 / 180 * pi, 45 / 180 * pi, 2 * λ)
    # dipole=convert(FitzgeraldDipole{Float64},HertzDipole([0.0,0.0,0.0], [1.0,0.0,0.0], complex(1.0)))
    dipoles = [
        HertzDipole([0.0, 0.0, 2.5 * λ], [1.0, 0.0, 0.0], complex(1.0))
        HertzDipole([0.0, 0.0, 2.5 * λ], [0.0, 1.0, 0.0], complex(1.0))
        HertzDipole([0.0, 0.0, 2.5 * λ], [0.0, 0.0, 1.0], complex(1.0))
        HertzDipole([0.0, 2.5 * λ, 0.0], [1.0, 0.0, 0.0], complex(1.0))
        HertzDipole([0.0, 2.5 * λ, 0.0], [0.0, 1.0, 0.0], complex(1.0))
        HertzDipole([0.0, 2.5 * λ, 0.0], [0.0, 0.0, 1.0], complex(1.0))
        HertzDipole([2.5 * λ, 0.0, 0.0], [1.0, 0.0, 0.0], complex(1.0))
        HertzDipole([2.5 * λ, 0.0, 0.0], [0.0, 1.0, 0.0], complex(1.0))
        HertzDipole([2.5 * λ, 0.0, 0.0], [0.0, 0.0, 1.0], complex(1.0))
        FitzgeraldDipole([0.0, 0.0, 2.5 * λ], [1.0, 0.0, 0.0], complex(1.0))
        FitzgeraldDipole([0.0, 0.0, 2.5 * λ], [0.0, 1.0, 0.0], complex(1.0))
        FitzgeraldDipole([0.0, 0.0, 2.5 * λ], [0.0, 0.0, 1.0], complex(1.0))
        FitzgeraldDipole([0.0, 2.5 * λ, 0.0], [1.0, 0.0, 0.0], complex(1.0))
        FitzgeraldDipole([0.0, 2.5 * λ, 0.0], [0.0, 1.0, 0.0], complex(1.0))
        FitzgeraldDipole([0.0, 2.5 * λ, 0.0], [0.0, 0.0, 1.0], complex(1.0))
        FitzgeraldDipole([2.5 * λ, 0.0, 0.0], [1.0, 0.0, 0.0], complex(1.0))
        FitzgeraldDipole([2.5 * λ, 0.0, 0.0], [0.0, 1.0, 0.0], complex(1.0))
        FitzgeraldDipole([2.5 * λ, 0.0, 0.0], [0.0, 0.0, 1.0], complex(1.0))
    ]
    for dipole in dipoles
        αHertz = convertrepresentation(IncidentSphericalExpansion{ComplexF64}, [dipole], k0)
        for obsdip in Hertzobsdips
            E1 = efield([dipole], obsdip.pos, k0)
            H1 = hfield([dipole], obsdip.pos, k0)
            E2, H2 = ehfield(αHertz, obsdip.pos, k0)
            @test all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
            @test all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))
            # println(E1)
            # println(E2)
            # println(E1./E2)

        end

        xvec = -1.25λ:0.75*λ:1.25λ
        yvec = -1.25λ:0.75*λ:1.25λ
        zvec = -1.25λ:0.75*λ:1.25λ

        for x in xvec, y in yvec, z in zvec
            E1 = efield([dipole], [x; y; z], k0)
            H1 = hfield([dipole], [x; y; z], k0)
            E2, H2 = ehfield(αHertz, [x; y; z], k0)
            @test all(isapprox.(E1, E2; rtol = 1e-4, atol = 1e-6))
            @test all(isapprox.(H1, H2; rtol = 1e-4, atol = 1e-6))

        end
    end
end

