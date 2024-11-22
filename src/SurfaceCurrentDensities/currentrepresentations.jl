
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
        deepcopy(s.functionspace),
        similar(s.excitations),
        s.wavenumber,
    )
end
function getwavenumber(currents::SurfaceCurrentDensity)
    return currents.wavenumber
end
function setwavenumber!(
    currents::SurfaceCurrentDensity{P,E,S,C},
    val::Number,
) where {P,E,S,C}
    return SurfaceCurrentDensity{P,E,S,C}(currents.functionspace, currents.excitations, val)
end
function Base.getindex(currents::SurfaceCurrentDensity, i)
    return getindex(currents.excitations, i)
end
function Base.setindex!(currents::SurfaceCurrentDensity, i, v)
    return setindex!(currents.excitations, i, v)
end
function Base.size(currents::SurfaceCurrentDensity)
    return size(currents.excitations)
end

function numfunctions(currentdensity::SurfaceCurrentDensity)
    return BEAST.numfunctions(functionspace(currentdensity))
end

function equivalentorder(currentdensity::SurfaceCurrentDensity; ϵ = 1e-7)
    L = _modeorder(rsph, k0; ϵ = ϵ)
    return L
end

#### Field evaluations
function farfield(
    currents::SurfaceCurrentDensity{P,E,B,C},
    θϕ::Tuple{T,T},
    k₀::Number,
) where {P<:PropagationType,E<:ElmagType,B<:BEAST.Space{<:Real},C,T}
    θ, ϕ = θϕ

    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)

    pts = [point(cosϕ * sinθ, sinϕ * sinθ, cosθ)]
    ffd = potential(
        MWFarField3D(; wavenumber = k₀),
        pts,
        currents.excitations,
        currents.functionspace,
    )[1]
    Eθ = convert.(Complex{T}, udot([cosθ * cosϕ, cosθ * sinϕ, -sinθ], ffd))
    Eϕ = convert.(Complex{T}, udot([-sinϕ, cosϕ, zero(eltype(θvec))], ffd))

    return _weightedfarfieldpolarization(E(), Eθ, Eϕ)
end

function _weightedfarfieldpolarization(::Electric, Eθ, Eϕ)
    factor = complex(0.0, -k₀) * Z₀ / (4π)
    return factor * Eθ, factor * Eϕ
end
function _weightedfarfieldpolarization(::Magnetic, Eθ, Eϕ)
    factor = complex(0.0, -k₀) / (4π)
    return factor * Eϕ, -factor * Eθ
end

function efield!(
    storage,
    currents::SurfaceCurrentDensity{Radiated,Electric,S,C},
    R::AbstractVector{C},
) where {C<:Number,S}
    gridpoint = [point(R[1], R[2], R[3])]
    store(v, m, n) = (storage .+= v * coeffs[n])
    potential!(
        store,
        BEAST.MWSingleLayerField3D(; wavenumber = getwavenumber(currents)),
        gridpoint,
        currents.excitations * Z₀,
        currents.functionspace,
    )
    return storage
end
function hfield!(
    storage,
    currents::SurfaceCurrentDensity{Radiated,Magnetic,S,C},
    R::AbstractVector{C},
) where {C<:Number,S}
    gridpoint = [point(R[1], R[2], R[3])]
    store(v, m, n) = (storage .+= v * coeffs[n])
    potential!(
        store,
        BEAST.MWDoubleLayerField3D(; wavenumber = getwavenumber(currents)),
        gridpoint,
        -currents.excitations,
        currents.functionspace,
    )
    return storage
end

#### Transformations
# function rotate(currents::SurfaceCurrentDensity, χ::Number, θ::Number, ϕ::Number)
#     newcurrents = deepcopy(currents)
#     if abs(χ) > 1e-16
#         CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, 0.0, χ])
#     end
#     if abs(θ) > 1e-16
#         CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, θ, 0.0])
#     end
#     if abs(ϕ) > 1e-16
#         CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, 0.0, ϕ])
#     end
#     return newcurrents
# end
function rotate!(
    rotated_currents::SurfaceCurrentDensity,
    currents::SurfaceCurrentDensity,
    χ::Number,
    θ::Number,
    ϕ::Number,
)

    # for i in eachindex(rotated_currents.functionspace.geo.vertices)
    #     rotated_currents.functionspace.geo.vertices[i] .= currents.functionspace.geo.vertices[i]
    # end
    rotated_currents.functionspace.geo.vertices .= currents.functionspace.geo.vertices
    rotated_currents.excitations .= currents.excitations



    if abs(χ) > 1e-16
        CompScienceMeshes.rotate!(rotated_currents.functionspace.geo, [0.0, 0.0, χ])
    end
    if abs(θ) > 1e-16
        CompScienceMeshes.rotate!(rotated_currents.functionspace.geo, [0.0, θ, 0.0])
    end
    if abs(ϕ) > 1e-16
        CompScienceMeshes.rotate!(rotated_currents.functionspace.geo, [0.0, 0.0, ϕ])
    end
    return rotated_currents
end

# function spatialshift(currents::SurfaceCurrentDensity, R)
#     newcurrents = deepcopy(currents)
#     CompScienceMeshes.translate!(newcurrents.functionspace.geo, -R)
#     return newcurrents
# end
function spatialshift!(
    shifted_currents::SurfaceCurrentDensity,
    currents::SurfaceCurrentDensity,
    R,
)
    shifted_currents.functionspace.geo.vertices .= currents.functionspace.geo.vertices
    shifted_currents.excitations .= currents.excitations
    CompScienceMeshes.translate!(shifted_currents.functionspace.geo, -R)
    return shifted_currents
end
