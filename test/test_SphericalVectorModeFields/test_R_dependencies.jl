using LinearAlgebra

import AntennaFieldRepresentations.R_dependencies_array
import AntennaFieldRepresentations.zc_ℓ
import AntennaFieldRepresentations.oneoverkA_deriv_zc_ℓ

Lmax=7

kAlist=[0.01, 3.57, 1000.1]
clist=[IncidentSphericalExpansion{ComplexF64}, AbsorbedSphericalExpansion{ComplexF64},RadiatingSphericalExpansion{ComplexF64}]
for kA in kAlist
    for k=0:Lmax
        if k< 2*kA+1
            
            # println(kA)
            for ck=1:3
                c=clist[ck]        
                zℓ, dzℓ= R_dependencies_array(c, Lmax+1, kA)
                
                global goal=zc_ℓ(c,k,kA)
                @test norm(zℓ[k+1]-goal)/norm(goal)< 1e-10
                # println(goal)
                # println(zn[k+1])
                # println("\n")
                goal=oneoverkA_deriv_zc_ℓ(c,k,kA)
                @test norm(dzℓ[k+1]-goal)/norm(goal)< 1e-10
                # println(goal)
                # println(dzn[k+1])
                # println("\n")
            # println("\n")
            end
        end
    end
end