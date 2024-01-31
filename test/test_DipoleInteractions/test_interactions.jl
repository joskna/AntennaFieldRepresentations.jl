using LinearAlgebra
# const c₀ = 299792458 # m/s
# const ε₀ = 8.854187812813e-12 # F/M
# const μ₀ = 1.2566370621219e-6 # N/A²
# const Z₀ = sqrt(μ₀ / ε₀) # Ω
# c0 = 299792458
        fvec = [1.5e9, 706.2e6]
        for kf in eachindex(fvec)
            f= fvec[kf]
            λ = c₀ / f
            k0 = 2 * pi / λ

            #check for symmetryf
            R1 = [5.0, 0.0, 0.0]
            R2 = [0.0, 5.0, 0.0]
            b = transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole(R1, [0.0, 0.0, 1.0], 1.0 + 0.0im), k0)
            b2 = transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole(R2, [0.0, 0.0, 1.0], 1.0 + 0.0im), k0)
            @test b ≈ b2
            #check for input types
            b2=transmission(HertzDipole(complex.([0.0, 0.0, 0.0]), [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole(R1, [0.0, 0.0, 1.0], 1.0 + 0.0im), k0)
            @test b ≈ b2
            b2=transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole(complex.(R1), [0.0, 0.0, 1.0], 1.0 + 0.0im), k0)
            @test b ≈ b2
            b2=transmission(HertzDipole(complex.([0.0, 0.0, 0.0]), [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole(complex.(R1), [0.0, 0.0, 1.0], 1.0 + 0.0im), k0)
            @test b ≈ b2

            sourcelist=[HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole([0.0, 0.0, 0.5], [0.0, 0.0, 1.0], 1.0 + 0.0im), FitzgeraldDipole([0.0, 0.0, 1.0], [0.0, 1.0, 0.0], 1.0 + 0.0im)]
            probelist=[HertzDipole([3*λ, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole([3*λ, 0.0, 0.5], [0.0, 0.0, 1.0], 1.0 + 0.0im), FitzgeraldDipole([3*λ, 0.0, 1.0], [0.0, 1.0, 0.0], 1.0 + 0.0im)]

            
            # test DipoleInteractionMatrix and its adjoint
            # import.AntennaDipoleInteractions.DipoleInteractionMatrix
            A=DipoleInteractionMatrix(sourcelist,probelist, k0)
            @test A[1:2,2:3]==[transmission(sourcelist[2],probelist[1],k0) transmission(sourcelist[3],probelist[1],k0);
            transmission(sourcelist[2],probelist[2],k0) transmission(sourcelist[3],probelist[2],k0) ]
                               
            v=[complex(1.0), complex(0.0, 1.0), complex(1.0,1.0)]
            b=A*v
            b2=zeros(ComplexF64, size(A,1))
            for i in eachindex(b2)
                for j in eachindex(v)
                    b2[i]+=transmission(sourcelist[j],probelist[i],k0)*v[j]
                end
            end
            b=adjoint(A)*v
            b2=zeros(ComplexF64, size(A,2))
            for i in eachindex(b2)
                for j in eachindex(v)
                    b2[i]+=conj.(transmission(sourcelist[i],probelist[j],k0))*v[j]
                end
            end
            @test b2==b

            
            #check against Hansen: "Spherical Near-Field Measurements" page 328
            θvec=(0:36:180)/180*π
            ϕvec=(0:36:359)/180*π
            rvec=[0.7*λ; 2*λ;13.3*λ;101.1*λ]
            # import AntennaDipoleInteractions.Z₀
            for kr in eachindex(rvec)
                r=rvec[kr]
                for kθ in eachindex(θvec)
                    st, ct=sincos(θvec[kθ])
                    for kϕ in eachindex(ϕvec)
                        sp, cp=sincos(ϕvec[kϕ])
                        eᵣ=[st*cp; st*sp; ct]
                        eθ=[ct*cp; ct*sp; -st]
                        eϕ=[-sp; +cp; 0.0]

                        pos=r*[st*cp; st*sp; ct]
                        exc=-√π*√6/(sqrt(Z₀)*k0) # unit excitation to radiate power of 1/2 W  

                        F201eᵣ=-√6/(2*√π)*cis(-k0*r)/(k0^2*r^2)*(1-1.0im/(k0*r))*ct
                        F201eᵣ*=k0*√Z₀
                        b=transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)), HertzDipole(pos, eᵣ, 1.0 + 0.0im), k0)
                        @test b≈0.5 * F201eᵣ

                        F201eθ=-√6/(4*√π)*cis(-k0*r)/(k0*r)*(1.0im+1/(k0*r)-1.0im/(k0^2*r^2))*st
                        F201eθ*=k0*√Z₀
                        b=transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)), HertzDipole(pos, eθ, 1.0 + 0.0im), k0)
                        @test b≈0.5 * F201eθ

                        F101eϕ=-√6/(4*√π)*cis(-k0*r)/(k0*r)*(1.0-1.0im/(k0*r))*st
                        F101eϕ*=1im*k0/(√Z₀)
                        b=transmission(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)), FitzgeraldDipole(pos, eϕ, -1.0 + 0.0im), k0)#test magnetic field with negative Fitzgerald dipole
                        @test b≈0.5 * F101eϕ atol=1e-15

                        F101eϕ=-√6/(4*√π)*cis(-k0*r)/(k0*r)*(1.0-1.0im/(k0*r))*st
                        F101eϕ*=1im*k0/(√Z₀)
                        b=transmission(FitzgeraldDipole(pos, eϕ, -1.0 + 0.0im),HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)), k0)
                        @test b≈0.5 * F101eϕ atol=1e-15
                        

                    end
                end
            end
        end

dipole=HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im)
fdipole=converttype(FitzgeraldDipole{eltype(dipole.pos)}, dipole)

@test fdipole.pos ≈ dipole.pos
@test fdipole.dir ≈ dipole.dir
@test fdipole.mag ≈ dipole.mag

@test typeof(fdipole)==FitzgeraldDipole{Float64, ComplexF64}

diplist=[HertzDipole([3.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole([3.0, 0.0, 0.5], [0.0, 0.0, 1.0], 1.0 + 0.0im), HertzDipole([3.0, 0.0, 1.0], [0.0, 1.0, 0.0], 1.0 + 0.0im)]
flist=converttype.(FitzgeraldDipole{eltype(dipole.pos)}, diplist)
@test typeof(flist)==Vector{FitzgeraldDipole{Float64, ComplexF64}}