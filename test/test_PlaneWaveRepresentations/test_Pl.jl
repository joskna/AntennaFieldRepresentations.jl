# Check against Wikipedia

import AntennaFieldRepresentations.collectPl
Lmax=10

xlist=[1e-4, 0.3,0.5,0.777, 1-1e-4, 1.00]
for k=1:length(xlist)
    x=xlist[k]
    Pℓ=collectPl(Lmax,x)

    goalval=1.0
    @test abs(Pℓ[1]-goalval)/abs(goalval)<1e-10

    goalval=x
    @test abs(Pℓ[2]-goalval)/abs(goalval)<1e-10 

    goalval=0.5*(3*x^2-1)
    @test abs(Pℓ[3]-goalval)/abs(goalval)<1e-10

    goalval=0.5*(5*x^3-3*x)
    @test abs(Pℓ[4]-goalval)/abs(goalval)<1e-10

    goalval=1/8*(35*x^4-30*x^2+3)
    @test abs(Pℓ[5]-goalval)/abs(goalval)<1e-10 

    goalval=1/8*(63*x^5-70*x^3+15*x)
    @test abs(Pℓ[6]-goalval)/abs(goalval)<1e-10

    goalval=1/16*(231*x^6-315*x^4+105*x^2-5)
    @test abs(Pℓ[7]-goalval)/abs(goalval)<1e-10

    goalval=1/16*(429*x^7-693*x^5+315*x^3-35*x)
    @test abs(Pℓ[8]-goalval)/abs(goalval)<1e-10

    goalval=1/128*(6435*x^8-12012*x^6+6930*x^4-1260*x^2+35)
    @test abs(Pℓ[9]-goalval)/abs(goalval)<1e-10

    goalval=1/128*(12155*x^9-25740*x^7+18018*x^5-4620*x^3+315*x)
    @test abs(Pℓ[10]-goalval)/abs(goalval)<1e-10

    goalval=1/256*(46189*x^10-109395*x^8+90090*x^6-30030*x^4+3465*x^2-63)
    @test abs(Pℓ[11]-goalval)/abs(goalval)<1e-10
end