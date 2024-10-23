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

@testset "Interaction Dipoles <-> Spherical" begin
    rad_aut = 1.51 * λ

    rad_meas = 1.5 * rad_aut

    x_aut = sqrt(rad_aut^2 - (λ / 4)^2) / sqrt(2)
    y_aut = x_aut

    Hertzobsdips = dipoles_on_sphere(12 / 180 * pi, 12 / 180 * pi, 2.5 * λ)
    Fitzobsdips = converttype.(FitzgeraldDipole{Float64}, Hertzobsdips)
    Hertzsrcdips = generate_AUTdips(
        collect(-λ/4:λ/4:λ/4),
        collect(-2x_aut/2:λ/4:2x_aut/2),
        collect(-x_aut/2:λ/4:x_aut/2),
        k0,
    )
    Fitzsrcdips = converttype.(FitzgeraldDipole{Float64}, Hertzsrcdips)


    bdips = complex(0.0)
    for k ∈ eachindex(Hertzobsdips), kk ∈ eachindex(Hertzsrcdips)
        bdips += transmission(Hertzsrcdips[kk], Hertzobsdips[k], k0)
    end
    bsph = transmission(
        Hertzsrcdips,
        convertrepresentation(IncidentSphericalExpansion{ComplexF64}, Hertzobsdips, k0),
        k0,
    )
    @test isapprox(bdips, bsph, rtol = 1e-8)
    bsph = transmission(
        convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, Hertzsrcdips, k0),
        Hertzobsdips,
        k0,
    )
    @test isapprox(bdips, bsph, rtol = 1e-8)

    bdips = complex(0.0)
    for k ∈ eachindex(Hertzobsdips), kk ∈ eachindex(Fitzsrcdips)
        bdips += transmission(Fitzsrcdips[kk], Hertzobsdips[k], k0)
    end
    bsph = transmission(
        Fitzsrcdips,
        convertrepresentation(IncidentSphericalExpansion{ComplexF64}, Hertzobsdips, k0),
        k0,
    )
    @test isapprox(bdips, bsph, rtol = 1e-8)

    bdips = complex(0.0)
    for k ∈ eachindex(Fitzobsdips), kk ∈ eachindex(Hertzsrcdips)
        bdips += transmission(Hertzsrcdips[kk], Fitzobsdips[k], k0)
    end
    bsph = transmission(
        convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, Hertzsrcdips, k0),
        Fitzobsdips,
        k0,
    )
    @test isapprox(bdips, bsph, rtol = 1e-8)

end


# @testset verbose=true "Interaction Dipole -> RadiatingSpherical" begin
#     rad_aut = 1.251* λ

# rad_meas =1.51*rad_aut

# x_aut = sqrt(rad_aut^2 - (λ / 4)^2) / sqrt(2)
# y_aut = x_aut



# autdips = generate_AUTdips( collect(-λ / 4:λ / 4:λ/4),
# collect(-2x_aut / 2:λ / 3:2x_aut / 2),
# collect(-x_aut / 2:λ / 3:x_aut / 2),
# k0)

# α=convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, autdips, k0)

# Hertzobsdips=dipoles_on_sphere(45/180*pi,45/180*pi,5*λ)
# Fitzobsdips=converttype.(FitzgeraldDipole{Float64},Hertzobsdips)
# for dipole in Hertzobsdips
#     b0=zero(ComplexF64)
#     for sourcedipole in autdips
#     b0+=transmission(sourcedipole,dipole,k0 )
#     end
#     b1=transmission(α,[dipole],k0)
#     b2=transmission([dipole],α,k0)
#     @test isapprox(b1 , b2; rtol=1e-6, atol=1e-6)
#     @test isapprox(b0 , b1; rtol=1e-6, atol=1e-6)
#     @test isapprox(b0 , b2; rtol=1e-6, atol=1e-6)

# end
# for dipole in Fitzobsdips
#     b0=zero(ComplexF64)
#     for sourcedipole in autdips
#     b0+=transmission(sourcedipole,dipole,k0 )
#     end
#     b1=transmission(α,[dipole],k0)
#     b2=transmission([dipole],α,k0)
#     @test isapprox(b1 , b2; rtol=1e-6, atol=1e-6)
#     @test isapprox(b0 , b1; rtol=1e-6, atol=1e-6)
#     @test isapprox(b0 , b2; rtol=1e-6, atol=1e-6)

# end


# end

