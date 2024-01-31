#Compare against Hansen p. 330
c0 = 299792458
f = 1.0e9

λ = c0 / f
k0 = 2 * pi / λ
Nmax=1
Jmax= Int((2 * ceil(Nmax) * (ceil(Nmax) + 2)))

ϑlist=[pi/6]
φlist=[0.0,pi/6, 0.9999*2*pi,2*pi ]

import AntennaFieldRepresentations.K_sℓm_array

for ϑ in ϑlist, φ in φlist
    global goal
# for ik= 1:length(ϑlist)
#     for ikk= 1:length(φlist)


#         global ϑ, goal
#         ϑ=ϑlist[ik]
#         φ=φlist[ikk]

        Kϑ, Kφ=K_sℓm_array(Jmax,ϑ,φ)

        j=sℓm_to_j(1,1,-1)
         goal=√3/(2)*exp(-1im*φ)/√(4π) 
        @test abs(Kϑ[j]-1im*goal)/abs(goal)<1e-10
        @test abs(Kφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

        j=sℓm_to_j(1,1,0)
         goal=-√6/(2)*sin(ϑ)/√(4π) 
        @test abs(Kϑ[j])<1e-10 # should be zero
        @test abs(Kφ[j]- goal)/abs(goal)<1e-10

        j=sℓm_to_j(1,1,1)
         goal=√3/(2)*exp(1im*φ)/√(4π) 
        @test abs(Kϑ[j]-1im*goal)/abs(goal)<1e-10
        @test abs(Kφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


        j=sℓm_to_j(2,1,-1)
        goal=-1im*√3/(2)*exp(-1im*φ)/√(4π) 
       @test abs(Kϑ[j]-(-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10
       @test abs(Kφ[j]- 1im*goal)/abs(goal)<1e-10

       j=sℓm_to_j(2,1,0)
        goal=-1im*√6/(2)*sin(ϑ)/√(4π) 
       @test abs(Kϑ[j]- goal)/abs(goal)<1e-10
       @test abs(Kφ[j])<1e-10 # should be zero

       j=sℓm_to_j(2,1,1)
        goal=-1im*√3/(2)*exp(1im*φ)/√(4π) 
       @test abs(Kϑ[j]-(cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10
       @test abs(Kφ[j]- 1im*goal)/abs(goal)<1e-10
    # end
end

ϑlist=[ 0.0, pi]
φlist=[0.0,pi/6, 0.9999*2*pi,2*pi ]


for ϑ in ϑlist, φ in φlist
    global goal

        Kϑ, Kφ=K_sℓm_array(Jmax,ϑ,φ)

        j=sℓm_to_j(1,1,-1)
         goal=√3/(2)*exp(-1im*φ)/√(4π) 
        @test abs(Kϑ[j]-1im*goal)/abs(goal)<1e-10
        @test abs(Kφ[j]- (cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10

        j=sℓm_to_j(1,1,0)
         goal=-√6/(2)*sin(ϑ)/√(4π) 
        @test abs(Kϑ[j])<1e-10 # should be zero
        @test abs(Kφ[j]- goal)<1e-10

        j=sℓm_to_j(1,1,1)
         goal=√3/(2)*exp(1im*φ)/√(4π) 
        @test abs(Kϑ[j]-1im*goal)/abs(goal)<1e-10
        @test abs(Kφ[j]- (-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10


        j=sℓm_to_j(2,1,-1)
        goal=-1im*√3/(2)*exp(-1im*φ)/√(4π) 
       @test abs(Kϑ[j]-(-cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10
       @test abs(Kφ[j]- 1im*goal)/abs(goal)<1e-10

       j=sℓm_to_j(2,1,0)
        goal=-1im*√6/(2)*sin(ϑ)/√(4π) 
       @test abs(Kϑ[j]- goal)<1e-10
       @test abs(Kφ[j])<1e-10 # should be zero

       j=sℓm_to_j(2,1,1)
        goal=-1im*√3/(2)*exp(1im*φ)/√(4π) 
       @test abs(Kϑ[j]-(cos(ϑ)*goal))/abs(cos(ϑ)*goal)<1e-10
       @test abs(Kφ[j]- 1im*goal)/abs(goal)<1e-10
    # end
end