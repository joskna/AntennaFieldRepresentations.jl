"""
    ProbeAntenna{A}

Wrapper around an `AntennaFieldRepresentation` to indicate that it is used as a probe.

# Type Parameters
- `A <: AntennaFieldRepresentation`
"""
struct ProbeAntenna{A<:AntennaFieldRepresentation}
    aut_field::A
    probesize::Real
end
function getwavenumber(pa::ProbeAntenna)
    return getwavenumber(pa.aut_field)
end
"""
    getprobesize(p::ProbeAntenna)

Return the radius of the smallest sphere fitting around the physical dimensions of the probe antenna.
"""
function getprobesize(p::ProbeAntenna)
    return p.probesize
end
"""
    getprobesize(p::ProbeAntenna)

Change the radius of the smallest sphere fitting around the physical dimensions of the probe antenna.
"""
function setprobesize!(p::ProbeAntenna, probesize::Real)
    p = ProbeAntenna(p.aut_field, probesize)
    return p
end

function transmit(aut_field, p::ProbeAntenna{D}) where {D<:DipoleArray{Radiated}}
    return transmit(aut_field, p, [0, 0, 0])
end
function transmit(
    aut_field, p::ProbeAntenna{D}, R
) where {D<:DipoleArray{Radiated,Electric}}
    dipoles = p.aut_field

    result = zero(eltype(dipoles))
    for k in eachindex(dipoles.positions)
        E = efield(aut_field, dipoles.positions[k] + R)
        result += 0.5 * udot(E, dipoles.orientations[k]) * dipoles.dipolemoments[k]
    end
    return result
end
function transmit(
    aut_field, p::ProbeAntenna{D}, R
) where {D<:DipoleArray{Radiated,Magnetic}}
    dipoles = p.aut_field

    result = zero(eltype(dipoles))
    for k in eachindex(dipoles.positions)
        H = hfield(aut_field, dipoles.positions[k] + R)
        result += -0.5 * udot(H, dipoles.orientations[k]) * dipoles.dipolemoments[k]
    end
    return result
end
