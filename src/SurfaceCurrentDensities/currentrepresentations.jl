
"""
    SurfaceCurrentDensity{P <: PropagationType, E <: ElmagType, S <: BEAST.Space{<: Real} , C}  <: AntennaFieldRepresentation{P, C}

Representation of an electromagnetic field via equivalent surface currents.

Behaves like an `AbstractVector{C}` with extra context.
The type parameter `E` defines if electric or magnetic surface current density.
"""
struct SurfaceCurrentDensity{P<:PropagationType,E<:ElmagType,S<:BEAST.Space{<:Real},C} <:
       AntennaFieldRepresentation{P,C}
    functionspace::S
    excitations::Vector{C}
    wavenumber::Number
end

function asvector(s::SurfaceCurrentDensity)
    return s.excitations
end
function Base.similar(s::SurfaceCurrentDensity{P,E,S,C}) where {P,E,S,C}
    return SurfaceCurrentDensity{P,E,S,C}(
        s.functionspace,
        similar(s.excitations),
        s.wavenumber,
    )
end
