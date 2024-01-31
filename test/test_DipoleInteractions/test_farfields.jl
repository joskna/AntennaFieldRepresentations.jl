# c0 = 299792458
        fvec = [1.5e9, 706.2e6]
        for kf in eachindex(fvec)
            f= fvec[kf]

            λ = c₀ / f
            k0 = 2 * pi / λ
            θvec=(0:36:180)/180*π
            ϕvec=(0:36:359)/180*π
            for kθ in eachindex(θvec)
                θ=θvec[kθ]
                for kϕ in eachindex(ϕvec)
                    ϕ=ϕvec[kϕ]
                    exc=-√π*√6/(sqrt(Z₀)*k0) # unit excitation to radiate power of 1/2 W
                    K201eθ=-1.0im*√6/2*sin(θ)*(√Z₀)/(sqrt(4*pi))
                    Eθ, Eϕ= farfield(HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)),θ,ϕ,k0)
                    @test Eθ≈K201eθ #atol=1.0e-14
                    @test Eϕ≈0 atol=1.0e-14
                    Eθ, Eϕ= farfield(FitzgeraldDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], exc* Z₀),θ,ϕ,k0)
                    K101eϕ=complex(0.0,1.0)*√6/2*sin(θ)*(√Z₀)/(sqrt(4*pi))
                    @test Eϕ≈K101eϕ atol=1.0e-14
                    @test Eθ≈0 atol=1.0e-14
                end
            end 
        end

        freq= 1.5e9

        λ = c₀ / freq
        k0 = 2 * pi / λ

θvec=collect(0:36:180)/180*π
ϕvec=collect(0:36:359)/180*π

#calculate far field  of Hertzian dipole at origin
exc=-√π*√6/(sqrt(Z₀)*k0)
Eθ, Eϕ= farfield([HertzDipole([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], complex(exc)),
                  FitzgeraldDipole([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], complex(exc))],
                 θvec,
                 ϕvec,
                 k0)