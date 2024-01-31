using LinearAlgebra



# copy from ElectromagneticDipoleUtils to not depend on this library
function generate_AUTdips(xvec::Array{Float64,1}, yvec::Array{Float64,1}, zvec::Array{Float64,1}, k0::Float64)

    nx = length(xvec)
    ny = length(yvec)
    nz = length(zvec)

    ycenter = (maximum(yvec) + minimum(yvec)) / 2
    ysize = maximum([abs(maximum(yvec) - ycenter),abs(minimum(yvec) - ycenter)])

    ndips = nx * ny * nz

    dipoles = Array{HertzDipole{Float64},1}(undef, ndips)# Array(HertzDipole{Float64},ndips)
 # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index =  (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
           # println(index)
                pos = [xvec[kkk],yvec[kk],zvec[k]]
           # println(pos)
                dipoles[index] = HertzDipole(pos, [0.0 + 0.0im,0.0 + 0.0im,1.0 + 0.0im], mag/ndips)

            end
        end
    end
    return dipoles

end
# copy from ElectromagneticDipoleUtils to not depend on this library
function dipoles_on_sphere(deltatheta::Float64, deltaphi::Float64, radius::T; tangential=true, pole=false) where T<:Number
    phi = 0:deltaphi:2 * pi - deltaphi;
    theta = deltatheta/2:deltatheta:pi - deltatheta/2
    if pole==true
        theta = 0:deltatheta:pi;
    end


    numtheta = length(theta)
    numphi = length(phi)

    fac=1/(numtheta*numphi)

    dipoles=Array{HertzDipole{typeof(radius)},1}(undef,2 * numtheta * numphi)
    # two tangential dipoles at each location

    if tangential==false
     # return also  radial dipoles  
    dipoles=Array{HertzDipole{typeof(radius)},1}(undef,3 * numtheta * numphi)
    end
    


    for k = 1:numtheta  
        cost=cos(theta[k])
        sint=sin(theta[k])

        for kk = 1:numphi
            cosp=cos(phi[kk])
            sinp=sin(phi[kk])

            index = numphi * (k - 1) + kk;
            x=radius * sint * cosp
            y=radius * sint * sinp
            z=radius*cost

            
            e_θ=[cost*cosp; cost*sinp; -sint]
            e_ϕ=[-sinp; cosp; 0.0]           
            dipoles[index]=HertzDipole([x;y;z], complex(e_θ), 1.0+0.0im)
            dipoles[index+numtheta * numphi ]=HertzDipole([x;y;z], complex(e_ϕ), fac+0.0im)
            if tangential== false
            e_r=[sint*cosp; sint*sinp; cost]
            dipoles[index+2*numtheta * numphi ]=HertzDipole([x;y;z], complex(e_r), fac+0.0im)
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


rad_aut = 1.51* λ

rad_meas =1.5*rad_aut

x_aut = sqrt(rad_aut^2 - (λ / 4)^2) / sqrt(2)
y_aut = x_aut
autdips=generate_AUTdips( collect(-λ / 4:λ / 4:λ/4),
collect(-2x_aut / 2:λ / 4:2x_aut / 2),
collect(-x_aut / 2:λ / 4:x_aut / 2),
k0)
autdips=rotate.(autdips; χ=pi/2.7, θ=pi/1.3, ϕ=-0.7*pi)
α=convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, autdips, k0)
_,Lmax,__=j_to_sℓm(length(α.coefficients))
FFdip=convertrepresentation(FarfieldPattern{ComplexF64},autdips, Lmax, k0)


@testset "Conversion Farfield <-> RadiatingSpherical" begin



ffsph=convertrepresentation(FarfieldPattern{ComplexF64},α,k0)
@test all(isapprox.(FFdip.Eθ,ffsph.Eθ; rtol=1e-8, atol=1e-7))
@test all(isapprox.(FFdip.Eϕ,ffsph.Eϕ; rtol=1e-8, atol=1e-7))

αret=convertrepresentation(RadiatingSphericalExpansion{ComplexF64}, ffsph)
@test norm(αret.coefficients- α.coefficients)/norm(α.coefficients) < 1e-9
@test all(isapprox.(αret.coefficients, α.coefficients; rtol=1e-8, atol=1e-8))
end

@testset "Conversion PlaneWave <-> IncidentSpherical" begin
    pw=PlaneWave(k0*[0.0,0.0,-1.0],[1.0,0.0,0.0],complex(1.0,0.0))
    αPW=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, pw,20)
    pws=convertrepresentation(PlaneWaveSpectrum{ComplexF64}, αPW, k0)
    αPWS=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, pws)
    αinc= AntennaFieldRepresentations.αinc_PW(20)/complex(0.0,1256.6370621200554)
    
    @test norm(αPW.coefficients- αPWS.coefficients)/norm(αPW.coefficients) < 1e-9
    @test all(isapprox.(αPWS.coefficients, αPW.coefficients; rtol=1e-8, atol=1e-8))

    @test norm(αPW.coefficients- αinc)/norm(αPW.coefficients) < 1e-9
    @test all(isapprox.(αinc, αPW.coefficients; rtol=1e-8, atol=1e-8))
end


@testset "Conversion PlanewaveSpectrum <-> IncidentSpherical" begin
R=[-10*λ,0.0,0.0]
autdipsshifted=deepcopy(autdips)
for dipole in autdipsshifted
    dipole.pos+=R
end

pwsinc=translate(FFdip, -R, k0)
αpwsinc=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, pwsinc)
# αincref=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, autdipsshifted, length(αpwsinc.coefficients), k0)
αincref=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, autdipsshifted, j_to_sℓm(length(αpwsinc.coefficients))[2], k0)

Jmax=sℓm_to_j(2,19,19)
@test norm(αpwsinc.coefficients[1:Jmax]- αincref.coefficients[1:Jmax])/norm(αincref.coefficients[1:Jmax]) < 1e-5

@test all(isapprox.(αpwsinc.coefficients[1:Jmax], αincref.coefficients[1:Jmax]; rtol=1e-5, atol=1e-5))

pwsconv=convertrepresentation(PlaneWaveSpectrum{ComplexF64}, αpwsinc, k0)
αconv=convertrepresentation(IncidentSphericalExpansion{ComplexF64}, pwsconv)

@test all(isapprox.(αconv.coefficients,αpwsinc.coefficients ))

ztest=-1.2λ:λ/5:1.2λ
ytest=-1.2λ:λ/5:1.2λ
xtest=-1.2λ:λ/5:1.2λ
for x in xtest, y in ytest, z in ztest
    Rvec=[x,y,z]
    E1, H1=ehfield(autdips, Rvec-R, k0)
    E2, H2=ehfield(pwsinc, Rvec, k0)
    E3, H3=ehfield(αincref, Rvec, k0)
    E4, H4=ehfield(αpwsinc, Rvec, k0)
    E5, H5=ehfield(pwsconv, Rvec, k0)
    # @test (all(isapprox.(E1,E2; rtol=1e-4, atol=1e-6)))
    # @test (all(isapprox.(H1,H2; rtol=1e-4, atol=1e-6)))

    @test (all(isapprox.(E1,E3; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(E1,E3; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(E1-E3)/norm(E1))
        println(norm(E1))
        println(E1)
        println(E3)
    end
    @test (all(isapprox.(H1,H3; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(H1,H3; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(H1-H3)/norm(H1))
        println(norm(H1))
        println(H1)
        println(H3)
    end

    @test (all(isapprox.(E1,E4; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(E1,E4; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(E1-E4)/norm(E1))
        println(norm(E1))
        println(E1)
        println(E4)
    end
    @test (all(isapprox.(H1,H4; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(H1,H4; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(H1-H4)/norm(H1))
        println(norm(H1))
        println(H1)
        println(H4)
    end

    @test (all(isapprox.(E1,E5; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(E1,E5; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(E1-E5)/norm(E1))
        println(norm(E1))
        println(E1)
        println(E5)
    end
    @test (all(isapprox.(H1,H5; rtol=1e-4, atol=1e-6)))
    if !all(isapprox.(H1,H5; rtol=1e-4, atol=1e-6))
        println(x)
        println(y)
        println(z)
        println(norm(H1-H5)/norm(H1))
        println(norm(H1))
        println(H1)
        println(H5)
    end

end
end

