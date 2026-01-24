
using LinearAlgebra

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = generate_AUTdips(
    collect((-0.25λ):(λ / 4):(0λ)),
    collect((-0.5λ):(λ / 4):(0.5λ)),
    collect((-1λ):(λ / 4):(1λ)),
    k0,
)
dipoles = AntennaFieldRepresentations.rotate(dipoles, 0.0, 0.9, 1.3)

pwe = changerepresentation(PlaneWaveExpansion, dipoles)

##########################################################################
# define measurement positions

# prototype probe
protopositions = [[0.0, 0.0, 0.0]]
magnitudes = [complex(1.0)]

####### spherical
probeθ = ProbeAntenna(
    HertzArray(
        protopositions,
        [complex.([1.0, 0.0, 0.0]) for k in 1:length(protopositions)],
        magnitudes,
        k0,
    ),
    λ / 10,
)
probeϕ = ProbeAntenna(
    HertzArray(
        protopositions,
        [complex.([0.0, 1.0, 0.0]) for k in 1:length(protopositions)],
        magnitudes,
        k0,
    ),
    λ / 10,
)
# measurement locations 
θs = (collect(0:10:180)) ./ 180 .* pi
ϕs = collect(0:10:359) ./ 180 .* pi
# ϕs = collect(0:15:90) ./ 180 .* pi

radius = 20λ

positions = [
    [
        radius * [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)] + [3 * λ, 0.0, 0.0] for θ in θs,
        ϕ in ϕs
    ]
    [
        radius * [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)] + [3 * λ, 0.0, 0.0] for θ in θs,
        ϕ in ϕs
    ]
]

orientationsco = [[[0.0, θ, ϕ] for θ in θs, ϕ in ϕs]; [[0.0, θ, ϕ] for θ in θs, ϕ in ϕs]]
# orientationscross = [[pi / 2, θ, ϕ] for θ in θs, ϕ in ϕs]
probeIDs = ones(Int64, size(positions))
probeIDs[(end ÷ 2 + 1):end, :] .*= 2
probes = [probeθ, probeϕ]
# probe = probeθ.aut_field

bref = Matrix{ComplexF64}(undef, size(positions))
for k in eachindex(positions)
    χ, θ, ϕ = orientationsco[k]
    tmp_probe = probes[probeIDs[k]]
    probe = AntennaFieldRepresentations.rotate(tmp_probe.aut_field, χ, θ, ϕ)
    bref[k] = transmit((dipoles), ProbeAntenna(probe, λ / 10), positions[k])
end
#######################

sampling = IrregularFieldSampling(positions, orientationsco, probeIDs, probes)
#########################################################

##########################################################
A = AntennaFieldRepresentations.MLFMMTransmitMap(
    dipoles, sampling, dipoles.wavenumber; mintranslationlevel=0, verbose=false
)

b = reshape(A * dipoles, size(bref))
@test maximum(abs.(b - bref)) ./ maximum(abs.(bref)) < 1e-4

#######################################################
using CompScienceMeshes
using BEAST
using IterativeSolvers

minrad = 1.4λ

# sphere_mesh = meshsphere(minrad, λ / 3)
filenamemesh = joinpath("testdata", "sphere_rad1p4.msh")
sphere_mesh = CompScienceMeshes.read_gmsh_mesh(filenamemesh)
Γ = BEAST.raviartthomas(sphere_mesh)
Γ_ = BEAST.buffachristiansen(sphere_mesh)

currents_el = SurfaceCurrentDensity{Radiated,Electric,typeof(Γ),ComplexF64}(
    Γ, ones(ComplexF64, numfunctions(Γ)), getwavenumber(dipoles)
)
currents_mag = SurfaceCurrentDensity{Radiated,Magnetic,typeof(Γ_),ComplexF64}(
    Γ_, ones(ComplexF64, numfunctions(Γ_)), getwavenumber(dipoles)
)

B = AntennaFieldRepresentations.MLFMMTransmitMap(
    currents_el, sampling, dipoles.wavenumber; verbose=false, expectedaccuracy=1e-3
)

C = AntennaFieldRepresentations.MLFMMTransmitMap(
    currents_mag, sampling, dipoles.wavenumber; verbose=false, expectedaccuracy=1e-3
)

D = [B C]
Dᴴ = adjoint(D)
DDᴴ = D * Dᴴ

bvec = A * dipoles
Dᴴb = Dᴴ * bvec

y = zeros(ComplexF64, size(bvec))
minres!(y, DDᴴ, bvec / norm(bvec); maxiter=50, verbose=true, abstol=1e-3)
x = Dᴴ * y * norm(bvec)

bvec2 = D * x
@test (norm(bvec2 - bvec) / norm(bvec)) < 1e-3

Bfarfield = ChangeRepresentationMap(
    typeof(pwe), B.sourcestruct; samplingstrategy=pwe.samplingstrategy
)
Cfarfield = ChangeRepresentationMap(
    typeof(pwe), C.sourcestruct; samplingstrategy=pwe.samplingstrategy
)

pwe2 = changerepresentation(PlaneWaveExpansion, dipoles)
pwe2 .*= 0

pwe3 = changerepresentation(PlaneWaveExpansion, dipoles)
pwe3 .*= 0

pwe2 .=
    Bfarfield * x[1:size(B, 2)] + Cfarfield * x[(size(B, 2) + 1):(size(B, 2) + size(C, 2))]

pwe3 .= pwe .- pwe2

@test (norm(pwe3) / norm(pwe)) < 1.1e-3
