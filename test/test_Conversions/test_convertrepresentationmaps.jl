using AntennaFieldRepresentations
using LinearAlgebra
function generate_AUTdips(
    xvec::Array{Float64,1},
    yvec::Array{Float64,1},
    zvec::Array{Float64,1},
    k0::Float64,
)

    nx = length(xvec)
    ny = length(yvec)
    nz = length(zvec)

    ycenter = (maximum(yvec) + minimum(yvec)) / 2
    ysize = maximum([abs(maximum(yvec) - ycenter), abs(minimum(yvec) - ycenter)])

    ndips = nx * ny * nz
    positions = Vector{Vector{Float64}}(undef, ndips)
    magnitudes = Vector{ComplexF64}(undef, ndips)

    # println(size(dipoles))
    for kkk = 1:nx
        dx = maximum(xvec) - xvec[kkk] # determine phase shift for radiation into x direction

        for kk = 1:ny
            dy = abs(yvec[kk] - ycenter) / ysize
            mag = complex(cos(dy) * exp(1im * dx * k0))
            for k = 1:nz


                index = (k - 1) * ny * nx + (kk - 1) * nx + kkk # dipole along z-axis
                # println(index)
                positions[index] = [xvec[kkk], yvec[kk], zvec[k]]
                magnitudes[index] = mag
                # println(pos)


            end
        end
    end
    return HertzArray(
        positions,
        [complex.([0.0, 0.0, 1.0]) for k = 1:length(positions)],
        magnitudes,
        k0,
    )

end

Z₀ = 376.730313669
f = 1.5e9
λ = AntennaFieldRepresentations.c₀ / f
k0 = 2 * pi / λ

dipoles = rotate(
    generate_AUTdips(
        collect(-0.5λ:λ/4:0.5λ),
        collect(-0.5λ:λ/4:0λ),
        collect(-1λ:λ/4:1λ),
        k0,
    ),
    0.7,
    0.9,
    1.3,
)

swe = changerepresentation(SphericalWaveExpansion{Radiated}, dipoles)
swe_ref = deepcopy(swe)
crm = ChangeRepresentationMap(PlaneWaveExpansion, swe)
pwe = crm.targetrepresentation
pwe .= crm * swe_ref
pwe_ref = deepcopy(pwe)
icrm = ChangeRepresentationMap(SphericalWaveExpansion, pwe)
swe_ = similar(swe)
swe_ .= icrm * pwe_ref

@test norm(swe_ .- swe_ref) / norm(swe_ref) < 2e-14

icrm_ = AntennaFieldRepresentations.inverse(crm)
swe_ .= icrm_ * pwe_ref

@test norm(swe_ .- swe_ref) / norm(swe_ref) < 2e-14

crm_ = AntennaFieldRepresentations.inverse(icrm)
pwe_ = similar(pwe)
pwe_ .= crm_ * swe_ref

@test norm(pwe_ - pwe_ref) / norm(pwe_ref) < 1e-16

