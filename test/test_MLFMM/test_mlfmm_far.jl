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
# dipoles= generate_AUTdips(collect(-0.5λ: λ/4: 0.5λ), collect(-0.5λ: λ/4: 0λ), collect(-1λ: λ/4: 1λ), k0)

stuff = MLFMMSource(dipoles, dipoles.wavenumber; verbose=false)

samplingtype = GaussLegendreθRegularϕSampling

θvec, ϕvec = samples(AntennaFieldRepresentations._standardsampling(samplingtype, 5))

θϕs = [(θvec[k], ϕvec[kk]) for k in eachindex(θvec), kk in eachindex(ϕvec)]

using LinearAlgebra

pwe = changerepresentation(PlaneWaveExpansion, dipoles)

AntennaFieldRepresentations._aggregate_to_farfield!(stuff, dipoles)

swe = changerepresentation(SphericalWaveExpansion, dipoles)

begin
    AntennaFieldRepresentations._aggregate_to_farfield!(stuff, dipoles)
    swe2 = changerepresentation(SphericalWaveExpansion, stuff.nodefarfields[1])
end

@test norm(swe .- swe2[1:length(swe)]) / norm(swe) < 3e-5

crm = ChangeRepresentationMap(typeof(pwe), stuff; samplingstrategy=pwe.samplingstrategy)
crm2 = ChangeRepresentationMap(SphericalWaveExpansion, pwe)

CRM = crm2 * crm

swe2 = CRM * dipoles
@test norm(swe .- swe2[1:length(swe)]) / norm(swe) < 3e-5

# Also adjoint and transpose operations should work: 

A = ChangeRepresentationMap(typeof(pwe), stuff; samplingstrategy=pwe.samplingstrategy)
Aᴴ = adjoint(A)

AAᴴ = A * Aᴴ
AᴴA = Aᴴ * A

using IterativeSolvers

b = Vector(pwe)

Aᴴb = Aᴴ * b

x = zeros(ComplexF64, size(b))

minres!(x, AAᴴ, b / norm(b); maxiter=100, verbose=true, abstol=1e-3)

@test norm(b - AAᴴ * x * norm(b)) / norm(b) < 1e-3
