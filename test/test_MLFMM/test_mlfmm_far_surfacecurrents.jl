using AntennaFieldRepresentations

# function generate_AUTdips(
#     xvec::Array{Float64,1},
#     yvec::Array{Float64,1},
#     zvec::Array{Float64,1},
#     k0::Float64,
# )

#     nx = length(xvec)
#     ny = length(yvec)
#     nz = length(zvec)

#     ycenter = (maximum(yvec) + minimum(yvec)) / 2
#     ysize = maximum([abs(maximum(yvec) - ycenter), abs(minimum(yvec) - ycenter)])

#     ndips = nx * ny * nz
#     positions = Vector{Vector{Float64}}(undef, ndips)
#     magnitudes = Vector{ComplexF64}(undef, ndips)

#     # println(size(dipoles))
#     for kkk = 1:nx
#         dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

#         for kk = 1:ny
#             dy = abs(yvec[kk] - ycenter) / ysize
#             mag = complex(cos(dy) * exp(1im * dx * k0))
#             for k = 1:nz

#                 index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
#                 # println(index)
#                 positions[index] = [xvec[kkk], yvec[kk], zvec[k]]
#                 magnitudes[index] = mag
#                 # println(pos)

#             end
#         end
#     end
#     return HertzArray(
#         positions,
#         [complex.([0.0, 0.0, 1.0]) for k = 1:length(positions)],
#         magnitudes,
#         k0,
#     )

# end

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = rotate(
    generate_AUTdips(
        collect((-0.5λ):(λ / 4):(0.5λ)),
        collect((-0.5λ):(λ / 4):(0λ)),
        collect((-1λ):(λ / 4):(1λ)),
        k0,
    ),
    0.7,
    0.9,
    1.3,
)

samplingtype = GaussLegendreθRegularϕSampling

θvec, ϕvec = samples(AntennaFieldRepresentations._standardsampling(samplingtype, 5))

θϕs = [(θvec[k], ϕvec[kk]) for k in eachindex(θvec), kk in eachindex(ϕvec)]

using LinearAlgebra

pwe = changerepresentation(PlaneWaveExpansion, dipoles)

using CompScienceMeshes
using BEAST

radius = 1.5 * λ

sphere_mesh = meshsphere(radius, λ / 5)
Γ = BEAST.raviartthomas(sphere_mesh)

currents = SurfaceCurrentDensity{Radiated,Electric,typeof(Γ),ComplexF64}(
    Γ, rand(ComplexF64, numfunctions(Γ)), dipoles.wavenumber
)
currentsmag = SurfaceCurrentDensity{Radiated,Magnetic,typeof(Γ),ComplexF64}(
    Γ, rand(ComplexF64, numfunctions(Γ)), dipoles.wavenumber
)

mlfmmsrc = MLFMMSource(currents, currents.wavenumber; verbose=false)
mlfmmsrcm = MLFMMSource(currentsmag, currents.wavenumber; verbose=false)

b = Vector(pwe)

A = ChangeRepresentationMap(typeof(pwe), mlfmmsrc; samplingstrategy=pwe.samplingstrategy)
Aᴴ = adjoint(A)

B = ChangeRepresentationMap(typeof(pwe), mlfmmsrcm; samplingstrategy=pwe.samplingstrategy)
Bᴴ = adjoint(B)

C = [A B]
# C= A

Cᴴ = adjoint(C)

CCᴴ = C * Cᴴ
CᴴC = Cᴴ * C

Cᴴb = Cᴴ * b

# y=zeros(ComplexF64,size(Cᴴb))
y = zeros(ComplexF64, size(b))
# y=CCᴴ* b / norm(b)

using IterativeSolvers
# gmres!(y, CᴴC , Cᴴb / norm(Cᴴb), verbose=true, restart=1, maxiter = 200,abstol=1e-6)
minres!(y, CCᴴ, b / norm(b); verbose=true, maxiter=100, abstol=1e-3)

# cg!(y, CCᴴ, b / norm(b), verbose=true, maxiter = 100, abstol=1e-3)

# x= y * norm(Cᴴb)
x = Cᴴ * y * norm(b)
pwe2 = similar(pwe)
pwe2 .= C * x

pwe3 = similar(pwe)
pwe3 .= pwe - pwe2

@test norm(pwe3) / norm(pwe) < 1e-3
