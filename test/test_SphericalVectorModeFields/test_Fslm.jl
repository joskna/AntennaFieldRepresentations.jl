#Compare against Hansen p. 328
c0 = 299792458
f = 1.0e9

λ = c0 / f
k0 = 2 * pi / λ

Nmax=1
Jmax= Int((2 * ceil(Nmax) * (ceil(Nmax) + 2)))

ϑlist=[0.0, pi, pi*0.7891]
φlist=[0.0,pi/6, 0.9999*2*pi,2*pi ]
rlist=[0.01, 2.7,10000]

import AntennaFieldRepresentations.F_sℓm_spherical_array

for ϑ in ϑlist, φ in φlist, r in rlist
    global goal
            kr=k0*r
            Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,Radiated(),r,ϑ,φ, k0)

            j=sℓm_to_j(1,1,-1)
            goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,0)
            goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j])<1e-10 #should be zero
            @test abs(Fφ[j]- (sin(ϑ)*goal))<1e-10 #should be zero

            j=sℓm_to_j(1,1,1)
            goal=√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


            j=sℓm_to_j(2,1,-1)
            goal=-√3/(2*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)<1e-10 #should be zero
            goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

            j=sℓm_to_j(2,1,0)
            goal=-√6/(2*√π) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*cos(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1im+1/kr-1im/(kr^2))*sin(ϑ)
            @test abs(Fϑ[j]-goal)<1e-10 #should be zero
            @test abs(Fφ[j])<1e-10 #should be zero

            j=sℓm_to_j(2,1,1)
            goal=√3/(2*√π)*exp(1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)<1e-10 #should be zero
            goal=-√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10


            Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,Incident(),r,ϑ,φ, k0)

            ckr=cos(kr)
            skr=sin(kr)
            j=sℓm_to_j(1,1,-1)
            goal=-√3/(4*√π)*exp(-1im*φ) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,0)
            goal=√6/(4*√π) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j])<1e-10 #should be zero
            @test abs(Fφ[j]- (sin(ϑ)*goal))<1e-10 #should be zero

            j=sℓm_to_j(1,1,1)
            goal=-√3/(4*√π)*exp(1im*φ) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


            j=sℓm_to_j(2,1,-1)
            goal=√3/(2*√π)*exp(-1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)<1e-10 #should be zero
            goal=√3/(4*√π)*exp(-1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

            j=sℓm_to_j(2,1,0)
            goal=√6/(2*√π) * 1/(kr^2)*(-ckr+skr/kr)*cos(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√6/(4*√π) * 1/kr*(skr+ckr/kr-skr/(kr^2))*sin(ϑ)
            @test abs(Fϑ[j]-goal)<1e-10 #should be zero
            @test abs(Fφ[j])<1e-10 #should be zero

            j=sℓm_to_j(2,1,1)
            goal=-√3/(2*√π)*exp(1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)<1e-10 #should be zero
            goal=-√3/(4*√π)*exp(1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10
end


# ϑlist=[ pi/6, pi/2, pi*0.7891]
# φlist=[0.0,pi/6, pi/2, 1.12345678*pi/2, 1.12345678*pi, 0.9999*2*pi,2*pi ]
# rlist=[0.01, 0.123, 2.7, 10, 100, 10000]

# for ϑ in ϑlist, φ in φlist, r in rlist
#     global goal
#             kr=k0*r
#             Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,RadiatingSphericalExpansion{ComplexF64},r,ϑ,φ, k0)

#             j=sℓm_to_j(1,1,-1)
#             goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
#             @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

#             j=sℓm_to_j(1,1,0)
#             goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1-1im/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j])<1e-10 #should be zero
#             @test abs(Fφ[j]- (sin(ϑ)*goal))/abs(sin(ϑ)*goal)<1e-10

#             j=sℓm_to_j(1,1,1)
#             goal=√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
#             @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


#             j=sℓm_to_j(2,1,-1)
#             goal=-√3/(2*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
#             @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
#             @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

#             j=sℓm_to_j(2,1,0)
#             goal=-√6/(2*√π) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*cos(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1im+1/kr-1im/(kr^2))*sin(ϑ)
#             @test abs(Fϑ[j]-goal)/abs(goal)<1e-10
#             @test abs(Fφ[j])<1e-10 #should be zero

#             j=sℓm_to_j(2,1,1)
#             goal=√3/(2*√π)*exp(1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=-√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
#             @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
#             @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10


#             Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,IncidentSphericalExpansion{ComplexF64},r,ϑ,φ, k0)

#             ckr=cos(kr)
#             skr=sin(kr)
#             j=sℓm_to_j(1,1,-1)
#             goal=-√3/(4*√π)*exp(-1im*φ) * 1/kr*(-ckr+skr/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
#             @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

#             j=sℓm_to_j(1,1,0)
#             goal=√6/(4*√π) * 1/kr*(-ckr+skr/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j])<1e-10 #should be zero
#             @test abs(Fφ[j]- (sin(ϑ)*goal))/abs(sin(ϑ)*goal)<1e-10

#             j=sℓm_to_j(1,1,1)
#             goal=-√3/(4*√π)*exp(1im*φ) * 1/kr*(-ckr+skr/kr)
#             @test abs(Fr[j])<1e-10 #should be zero
#             @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
#             @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


#             j=sℓm_to_j(2,1,-1)
#             goal=√3/(2*√π)*exp(-1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=√3/(4*√π)*exp(-1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
#             @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
#             @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

#             j=sℓm_to_j(2,1,0)
#             goal=√6/(2*√π) * 1/(kr^2)*(-ckr+skr/kr)*cos(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=-√6/(4*√π) * 1/kr*(skr+ckr/kr-skr/(kr^2))*sin(ϑ)
#             @test abs(Fϑ[j]-goal)/abs(goal)<1e-10
#             @test abs(Fφ[j])<1e-10 #should be zero

#             j=sℓm_to_j(2,1,1)
#             goal=-√3/(2*√π)*exp(1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
#             @test abs(Fr[j]-goal)/abs(goal)<1e-10
#             goal=-√3/(4*√π)*exp(1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
#             @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
#             @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10
# end

ϑlist=[ big(1e-4), big(pi-1e-4)]
φlist=[pi/6]
# φlist=[0.0,pi/6, pi/2, 2*pi ]
rlist=[2.7]

for ϑ in ϑlist, φ in φlist, r in rlist
    global goal
            kr=k0*r
            Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,Radiated(),r,ϑ,φ, k0)

            j=sℓm_to_j(1,1,-1)
            goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,0)
            goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j])<1e-10 #should be zero
            @test abs(Fφ[j]- (sin(ϑ)*goal))/abs(sin(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,1)
            goal=√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/kr*(1-1im/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


            j=sℓm_to_j(2,1,-1)
            goal=-√3/(2*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=√3/(4*√π)*exp(-1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

            j=sℓm_to_j(2,1,0)
            goal=-√6/(2*√π) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*cos(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√6/(4*√π) * exp(-1im*kr)/kr*(1im+1/kr-1im/(kr^2))*sin(ϑ)
            @test abs(Fϑ[j]-goal)/abs(goal)<1e-10
            @test abs(Fφ[j])<1e-10 #should be zero

            j=sℓm_to_j(2,1,1)
            goal=√3/(2*√π)*exp(1im*φ) * exp(-1im*kr)/(kr^2)*(1-1im/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√3/(4*√π)*exp(1im*φ) * exp(-1im*kr)/(kr)*(1im+1/kr-1im/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10


            Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,Incident(),r,ϑ,φ, k0)

            ckr=cos(kr)
            skr=sin(kr)
            j=sℓm_to_j(1,1,-1)
            goal=-√3/(4*√π)*exp(-1im*φ) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,0)
            goal=√6/(4*√π) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j])<1e-10 #should be zero
            @test abs(Fφ[j]- (sin(ϑ)*goal))/abs(sin(ϑ)*goal)<1e-10

            j=sℓm_to_j(1,1,1)
            goal=-√3/(4*√π)*exp(1im*φ) * 1/kr*(-ckr+skr/kr)
            @test abs(Fr[j])<1e-10 #should be zero
            @test abs(Fϑ[j]-1im*goal)/abs(goal)<1e-10
            @test abs(Fφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


            j=sℓm_to_j(2,1,-1)
            goal=√3/(2*√π)*exp(-1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=√3/(4*√π)*exp(-1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (-1im*goal))/abs(goal)<1e-10

            j=sℓm_to_j(2,1,0)
            goal=√6/(2*√π) * 1/(kr^2)*(-ckr+skr/kr)*cos(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√6/(4*√π) * 1/kr*(skr+ckr/kr-skr/(kr^2))*sin(ϑ)
            @test abs(Fϑ[j]-goal)/abs(goal)<1e-10
            @test abs(Fφ[j])<1e-10 #should be zero

            j=sℓm_to_j(2,1,1)
            goal=-√3/(2*√π)*exp(1im*φ) * 1/(kr^2)*(-ckr+skr/kr)*sin(ϑ)
            @test abs(Fr[j]-goal)/abs(goal)<1e-10
            goal=-√3/(4*√π)*exp(1im*φ) * 1/(kr)*(skr+ckr/kr-skr/(kr^2))
            @test abs(Fϑ[j]-cos(ϑ)*goal)/abs(cos(ϑ)*goal)<1e-10
            @test abs(Fφ[j]- (1im*goal))/abs(goal)<1e-10
end

Nmax=10
Jmax= Int((2 * ceil(Nmax) * (ceil(Nmax) + 2)))

ϑlist=[0.0im, pi/2+0.0im, pi/2+5.0im]
φlist=[0.0, pi/6]
# φlist=[0.0,pi/6, pi/2, 2*pi ]
rlist=[2.7]
for ϑ in ϑlist, φ in φlist, r in rlist
    global goal
            kr=k0*r
            Fr, Fϑ, Fφ=F_sℓm_spherical_array(Jmax,Radiated(),r,ϑ,φ, k0)
            for ikkkk=1:Jmax
                @test abs(Fr[ikkkk])<Inf
                @test abs(Fϑ[ikkkk])<Inf
                @test abs(Fφ[ikkkk])<Inf
            end

end