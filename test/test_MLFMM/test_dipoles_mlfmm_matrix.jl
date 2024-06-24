using AntennaFieldRepresentations
using BEAST, CompScienceMeshes
using StaticArrays
using ClusterTrees
using LinearAlgebra
# using PlotlyJS
using CompScienceMeshes

#### function definitions


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

# define physical constants
c = 3e8
μ0 = 4 * π * 1e-7
μr = 1.0
μ = μ0 * μr
ε0 = 8.854187812e-12
εr = 5
ε = ε0 * εr
c = 1 / sqrt(ε * μ)
f = 5e7
λ = c / f
κ = 2 * π / λ
ω = κ * c

# set MLFMM parameters
ϵ = 1e-4
minhalfsize=0.1*λ

#### Setup source

dmax=2.75*λ
xvec=collect(range(-dmax, dmax,24)) 
yvec=collect(range(-dmax/2, dmax/2,12)) 
zvec=collect(range(-λ/4, 0,2)) 
autdips= generate_AUTdips(zvec, xvec, yvec,  κ)
autdips_ref=deepcopy(autdips)
autdips_magone=deepcopy(autdips)
pos=Vector{SVector{3,Float64}}(undef, length(autdips))
for (index, dipole) in enumerate(autdips_ref)
    pos[index] = dipole.pos
    autdips_ref[index] = HertzDipole([0.0;0.0;0.0], dipole.dir, complex(1.0))
    autdips_magone[index] = HertzDipole(dipole.pos, dipole.dir, complex(1.0))
end

functionspacesrc=(sourcefunctions=autdips_ref, points=pos)

@time begin
    sourcestruct=AntennaFieldRepresentations.MLFMMSource(
    functionspacesrc,
    κ;
    verbose=false,
    minhalfsize=minhalfsize,
    θstencilsize=7,
    ϕstencilsize=7
    )
end
nothing

# zvec_obs=[1.5*dmax]
# obsdips= generate_AUTdips(xvec, yvec, zvec_obs, κ)
# for (index, dipole) in enumerate(obsdips)
#     obsdips[index].mag = 1
# end
obsdips=dipoles_on_sphere(30*pi/180, 30*pi/180, 1.5*dmax)

obsdips_ref=deepcopy(obsdips)
pos=Vector{SVector{3,Float64}}(undef, length(obsdips))
for (index, dipole) in enumerate(obsdips_ref)
    pos[index] = dipole.pos
    obsdips_ref[index] = HertzDipole([0.0;0.0;0.0], dipole.dir, dipole.mag)
end

functionspacerec=(sourcefunctions=obsdips_ref, points=pos)

@time begin
    receivestruct=AntennaFieldRepresentations.MLFMMReceive(
        functionspacerec,
        sourcestruct,
        verbose=false,
        θstencilsize=7,
        ϕstencilsize=7,
        num_bufferboxes=3
        )
end
nothing

x=[dipole.mag for dipole in autdips]
dipolematrix=DipoleInteractionMatrix(autdips_magone, obsdips, κ)
@time begin
    dipolematrix=DipoleInteractionMatrix(autdips_magone, obsdips, κ)
    b=dipolematrix*x
end

A= MLFMMInteractionMatrix{ComplexF64}(sourcestruct,receivestruct)
AᴴA=adjoint(A)*A
AAᴴ=A*adjoint(A)

@time begin
B=Base.materialize(dipolematrix)
C= Matrix{ComplexF64}(undef, size(A))
for i in 1:size(A,2)
    C[:,i] = A[:,i]
end
end

@time begin
D= Matrix{ComplexF64}(undef, size(adjoint(A)))
for i in 1:size(A,1)
    D[:,i] = adjoint(A)[:,i]
end
end

dev_forward=maximum(abs.(B-C)) / maximum(abs.(B))
dev_adjoint=maximum(abs.(adjoint(B)-D)) / maximum(abs.(B))

@test dev_forward < ϵ
@test dev_adjoint < ϵ