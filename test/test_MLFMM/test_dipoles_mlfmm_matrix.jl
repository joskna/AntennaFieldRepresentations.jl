
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
dipoles = rotate(dipoles, 0.0, 0.9, 1.3)

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
θs = (collect(0:35:180)) ./ 180 .* pi
ϕs = collect(0:55:359) ./ 180 .* pi
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
    probe = rotate(tmp_probe.aut_field, χ, θ, ϕ)
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
@test (maximum(abs.(b - bref)) ./ maximum(abs.(bref))) < 1.6e-5

Aᵀ = transpose(A)
Aᴴ = adjoint(A)

Amat = Matrix(A)
Amatᴴ = Matrix(Aᴴ)
Amatᵀ = Matrix(Aᵀ)

@test (maximum(abs.(Amat .- adjoint(Amatᴴ)))) < 1e-12
@test (maximum(abs.(Amat .- transpose(Amatᵀ)))) < 1e-12
