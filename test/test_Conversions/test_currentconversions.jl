using AntennaFieldRepresentations
using CompScienceMeshes
using BEAST
using LinearAlgebra
########
# Setup
#   |
#   v
########
vertices = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 1.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [0.5, 0.5, 1.0],
    [0.5, 1.0, 0.5],
    [1.0, 0.5, 0.5],
]

faces = [
    # front
    [1, 5, 10],
    [5, 6, 10],
    [6, 2, 10],
    [2, 1, 10],

    # down
    [1, 9, 5],
    [5, 9, 7],
    [7, 9, 3],
    [3, 9, 1],

    # left
    [1, 2, 11],
    [2, 4, 11],
    [4, 3, 11],
    [3, 1, 11],

    # right
    [5, 7, 14],
    [7, 8, 14],
    [8, 6, 14],
    [6, 5, 14],

    # back 
    [7, 13, 8],
    [8, 13, 4],
    [4, 13, 3],
    [3, 13, 7],

    # top
    [2, 6, 12],
    [6, 8, 12],
    [8, 4, 12],
    [4, 2, 12],
]

λ = 2.0
k0 = λ / (2π)

using StaticArrays
mesh = CompScienceMeshes.Mesh(SVector{3}.(vertices), SVector{3}.(faces))
Γ = BEAST.raviartthomas(mesh)

excitations = [
    -0.13 - 0.0014im
    0 - 1.6im
    0.13 + 1.6im
    -0.25 - 3.2im
    0 + 0im
    0.25 + 0
    0 + 1.6im
    -0.13 + 0.0014im
    0 + 0im
    -0.13 + 1.6im
    -0.25 + 3.2im
    0 - 1.6im
    -0.13 + 1.6im
    0 + 0im
    0.13 - 0.0014im
    -0.25 + 0im
    0.13 + 0.0014im
    0 - 1.6im
    0.13 + 1.6im
    0 + 0im
    0.13 + 0.0014im
    0.13 + 1.6im
    0 + 1.6im
    -0.25 + 0im
    -0.13 + 0.0014im
    -0.13 + 1.6im
    -0 + 1.6im
    -0.25 + 3.2im
    0.13 - 1.6im
    -0.13 + 0.0014im
    0 - 1.6im
    -0.25 + 0im
    0.13 + 0.0014im
    -0.13 - 1.6im
    -0.25 - 3.2im
    0 + 1.6im
]

currents = SurfaceCurrentDensity{Radiated,Electric,typeof(Γ),ComplexF64}(Γ, excitations, k0)

eldips = HertzArray([[0.5, 0.5, 0.5]], [complex.([0.0, 0.0, 1.0])], [complex(1.0)], k0)
magdips = FitzgeraldArray(
    [[0.5, 0.5, 0.5]],
    [complex.([0.0, -1.0, 0.0])],
    AntennaFieldRepresentations.Z₀ * [complex(1.0)],
    k0,
)
########
#   ^
#   |
# Setup
########

# convert into plane wave representation
pws_el = changerepresentation(PlaneWaveExpansion, eldips)
pws_mag = changerepresentation(PlaneWaveExpansion, magdips)

pws_dipoles = copy(pws_el)

pws_dipoles .= pws_el + pws_mag

pwe_currents = changerepresentation(
    PlaneWaveExpansion, currents; samplingstrategy=pws_dipoles.samplingstrategy
)

@test norm(pwe_currents - pws_dipoles) / norm(pws_dipoles) < 0.02

# convert into spherical wave representation
swe_currents = changerepresentation(SphericalWaveExpansion, currents)
swe_currents_Map = ChangeRepresentationMap(SphericalWaveExpansion, currents)

swe_currents2 = copy(swe_currents)
swe_currents2 .= 0

swe_currents2 .= swe_currents_Map * currents

@test norm(swe_currents2 - swe_currents) / norm(swe_currents) < 1e-15

swe_dipoles = changerepresentation(
    SphericalWaveExpansion, pws_dipoles; L=equivalentorder(swe_currents)
)
@test norm(swe_currents - swe_dipoles) ./ norm(swe_dipoles) < 0.02
