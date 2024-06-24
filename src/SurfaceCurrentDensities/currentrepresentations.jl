"""
Abstract type container type for Surface current densities
"""
abstract type AbstractSurfaceCurrentDensity <: AntennaFieldRepresentation end


struct ElectricSurfaceCurrentDensity{T<:Real, S<:BEAST.Space{T}} <: AbstractSurfaceCurrentDensity
    functionspace::S
    excitations::Vector{Complex{T}}
end
function ElectricSurfaceCurrentDensity(functionspace::BEAST.Space{T}) where{T}
    excitations=Vector{Complex{T}}(undef, length(functionspace.fns))
    return ElectricSurfaceCurrentDensity(functionspace, excitations)
end

struct MagneticSurfaceCurrentDensity{T<:Real, S<:BEAST.Space{T}} <: AbstractSurfaceCurrentDensity
    functionspace::S
    excitations::Vector{Complex{T}}
end
function MagneticSurfaceCurrentDensity(functionspace::BEAST.Space{T}) where{T}
    excitations=Vector{Complex{T}}(undef, length(functionspace.fns))
    return MagneticSurfaceCurrentDensity(functionspace, excitations)
end

struct ElmagSurfaceCurrentDensity{T<:Real, S<:BEAST.Space{T}} <: AbstractSurfaceCurrentDensity
    functionspace::S
    electricexcitations::Vector{Complex{T}}
    magneticexcitations::Vector{Complex{T}}
end
function ElmagSurfaceCurrentDensity(functionspace::BEAST.Space{T}) where{T}
    electricexcitations=Vector{Complex{T}}(undef, length(functionspace.fns))
    magneticexcitations=Vector{Complex{T}}(undef, length(functionspace.fns))
    return ElmagSurfaceCurrentDensity(functionspace, electricexcitations, magneticexcitations)
end

function functionspace(currentdensity)
    return currentdensity.functionspace
end
function numfunctions(currentdensity::AbstractSurfaceCurrentDensity)
    return BEAST.numfunctions(functionspace(currentdensity))
end
function numfunctions(currentdensity::ElmagSurfaceCurrentDensity)
    return 2 * BEAST.numfunctions(functionspace(currentdensity))
end

function elementtype(type::Type{ElectricSurfaceCurrentDensity{R, S}}) where {R<:Real,S}
    return Complex{R}
end
function elementtype(type::Type{MagneticSurfaceCurrentDensity{R,S}}) where {R<:Real,S}
    return Complex{R}
end
function elementtype(type::Type{ElmagSurfaceCurrentDensity{R,S}}) where {R<:Real,S}
    return Complex{R}
end

#### Type Conversions
function converttype(T::Type{ElectricSurfaceCurrentDensity{R,S}}, magneticcurrents::MagneticSurfaceCurrentDensity{R,S}) where {R<:Real,S}
    return ElectricSurfaceCurrentDensity(magneticcurrents.functionspace, magneticcurrents.excitations)
end

function converttype(T::Type{MagneticSurfaceCurrentDensity{R,S}}, electriccurrents::ElectricSurfaceCurrentDensity{R,S}) where {R<:Real,S}
    return MagneticSurfaceCurrentDensity(electriccurrents.functionspace, electriccurrents.excitations)
end




#### Field evaluations
function farfield(
    currents::ElectricSurfaceCurrentDensity{T,S}, θvec::Vector{<:Number}, ϕvec::Vector{<:Number}, k₀::Number
) where {T<:Real,S}

    nθ = length(θvec)
    nϕ = length(ϕvec)

    pts = [point(cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)) for ϕ in ϕvec for θ in θvec]
    ffd = reshape(potential(MWFarField3D(; wavenumber=k₀), pts, currents.excitations, currents.functionspace), (nθ, nϕ))

    Eθ =
        convert.(
            Complex{T},
            [
                udot([cos(θvec[kθ]) * cos(ϕvec[kϕ]), cos(θvec[kθ]) * sin(ϕvec[kϕ]), -sin(θvec[kθ])], ffd[kθ, kϕ]) for
                kθ in eachindex(θvec), kϕ in eachindex(ϕvec)
            ],
        )
    Eϕ =
        convert.(
            Complex{T},
            [udot([-sin(ϕvec[kϕ]), cos(ϕvec[kϕ]), zero(eltype(θvec))], ffd[kθ, kϕ]) for kθ in eachindex(θvec), kϕ in eachindex(ϕvec)],
        )

    factor = complex(zero(T), -k₀) * Z₀ / (4π)

    return factor * Eθ, factor * Eϕ
end

function farfield(
    currents::MagneticSurfaceCurrentDensity{T,S}, θvec::Vector{<:Number}, ϕvec::Vector{<:Number}, k₀::Number
) where {T<:Real,S}
    factor = 1 / Z₀
    Eθ, Eϕ = farfield(converttype(ElectricSurfaceCurrentDensity{T}, currents), θvec, ϕvec, k₀)
    return factor * Eϕ, -factor * Eθ
end

function farfield(currents::ElmagSurfaceCurrentDensity{T,S}, θ, ϕ, k₀::Number) where{T<:Real, S}
    Eθ1, Eϕ1 = farfield(ElectricSurfaceCurrentDensity(currents.functionspace, currents.electriceexcitations), θ, ϕ, k₀)
    Eθ2, Eϕ2 = farfield(MagneticSurfaceCurrentDensity(currents.functionspace, currents.magneticexcitations), θ, ϕ, k₀)
    return Eθ1 + Eθ2, Eϕ1 + Eϕ2
end


function efield(currents::ElectricSurfaceCurrentDensity{T,S}, Rvec::Array{<:AbstractVector{C}}, k₀::Number) where {T<:Real,C<:Number, S}
    gridpoints = [point(R[1], R[2], R[3]) for R in Rvec]
    E = potential(MWSingleLayerField3D(; wavenumber=k₀), gridpoints, currents.excitations * Z₀, currents.functionspace)
    return convert.(Array{Complex{T},1}, E)
end

function efield(currents::MagneticSurfaceCurrentDensity{T,S}, Rvec::Array{<:AbstractVector{C}}, k₀::Number) where {T<:Real,C<:Number,S}
    gridpoints = [point(R[1], R[2], R[3]) for R in Rvec]
    E = potential(BEAST.MWDoubleLayerField3D(; wavenumber=k₀), gridpoints, -currents.excitations, currents.functionspace)
    return convert.(Array{Complex{T},1}, E)
end


function hfield(currents::ElectricSurfaceCurrentDensity{T,S}, Rvec::Array{<:AbstractVector{C}}, k₀::Number) where {T<:Real,C<:Number,S}
    gridpoints = [point(R[1], R[2], R[3]) for R in Rvec]
    H = potential(BEAST.MWDoubleLayerField3D(; wavenumber=k₀), gridpoints, currents.excitations, currents.functionspace)
    return convert.(Array{Complex{T},1}, H)
end

function hfield(currents::MagneticSurfaceCurrentDensity{T,S}, Rvec::Array{<:AbstractVector{C}}, k₀::Number) where {T<:Real,C<:Number,S}
    gridpoints = [point(R[1], R[2], R[3]) for R in Rvec]
    H = potential(MWSingleLayerField3D(; wavenumber=k₀), gridpoints, currents.excitations * (1 / Z₀), currents.functionspace)
    return convert.(Array{Complex{T},1}, H)
end

function ehfield(currents::AbstractSurfaceCurrentDensity, R)
    return efield(currents, R), hfield(currents, R)
end

#### transmission
function transmission(
    sourcecurrents::ElectricSurfaceCurrentDensity{T,S1}, receivecurrents::ElectricSurfaceCurrentDensity{F,S2}, k₀::Number
) where {T<:Real,F<:Real,S1,S2}
    factor = complex(convert(T, 0.5 * Z₀))
    A = assemble(Maxwell3D.singlelayer(; wavenumber=k₀), sourcecurrents.functionspace, receivecurrents.functionspace)
    return factor * transpose(receivecurrents.excitations) * A * sourcecurrents.excitations
end

function transmission(
    sourcecurrents::ElectricSurfaceCurrentDensity{T,S1}, receivecurrents::MagneticSurfaceCurrentDensity{F,S2}, k₀::Number
) where {T<:Real,F<:Real,S1,S2}
    factor = -T(0.5)
    A = assemble(Maxwell3D.doublelayer(; wavenumber=k₀), sourcecurrents.functionspace, receivecurrents.functionspace)
    return factor * transpose(receivecurrents.excitations) * A * sourcecurrents.excitations
end

function transmission(
    sourcecurrents::MagneticSurfaceCurrentDensity{T,S1}, receivecurrents::ElectricSurfaceCurrentDensity{F,S2}, k₀::Number
) where {T<:Real,F<:Real,S1,S2}
    factor = -T(0.5)
    A = assemble(Maxwell3D.doublelayer(; wavenumber=k₀), sourcecurrents.functionspace, receivecurrents.functionspace)
    return factor * transpose(receivecurrents.excitations) * A * sourcecurrents.excitations
end

function transmission(
    sourcecurrents::MagneticSurfaceCurrentDensity{T,S1}, receivecurrents::MagneticSurfaceCurrentDensity{F,S2}, k₀::Number
) where {T<:Real,F<:Real,S1,S2}
    factor = complex(convert(T, 0.5 / Z₀))
    A = assemble(Maxwell3D.singlelayer(; wavenumber=k₀), sourcecurrents.functionspace, receivecurrents.functionspace)
    return factor * transpose(receivecurrents.excitations) * A * sourcecurrents.excitations
end

function transmission(
    sourcecurrents::AbstractSurfaceCurrentDensity, receivecurrents::ElmagSurfaceCurrentDensity{F,S}, k₀::Number
) where {F<:Real,S}
    elreccur = ElectricSurfaceCurrentDensity(receivecurrents.functionspace, receivecurrents.electricexcitations)
    magreccur = MagneticSurfaceCurrentDensity(receivecurrents.functionspace, receivecurrents.magneticexcitations)
    return transmission(sourcecurrents, elreccur, k₀) + transmission(sourcecurrents, magreccur, k₀)
end

function transmission(
    sourcecurrents::ElmagSurfaceCurrentDensity{F,S}, receivecurrents::AbstractSurfaceCurrentDensity, k₀::Number
) where {F<:Real,S}
    elsrccur = ElectricSurfaceCurrentDensity(sourcecurrents.functionspace, sourcecurrents.electricexcitations)
    magsrccur = MagneticSurfaceCurrentDensity(sourcecurrents.functionspace, sourcecurrents.magneticexcitations)
    return transmission(elsrccur, receivecurrents, k₀) + transmission(magsrccur, receivecurrents, k₀)
end

function transmission(
    sourcecurrents::ElmagSurfaceCurrentDensity{T,S1}, receivecurrents::ElmagSurfaceCurrentDensity{F,S2}, k₀::Number
) where {T<:Real,F<:Real,S1,S2}
    elsrccur = ElectricSurfaceCurrentDensity(sourcecurrents.functionspace, sourcecurrents.electricexcitations)
    magsrccur = MagneticSurfaceCurrentDensity(sourcecurrents.functionspace, sourcecurrents.magneticexcitations)
    return transmission(elsrccur, receivecurrents, k₀) + transmission(magsrccur, receivecurrents, k₀)
end

function transmission(field, receivecurrents::ElectricSurfaceCurrentDensity{R,S}, k₀::Number) where {R<:Real,S}
    incident_field = BEAST.assemble(((BEAST.n × (x -> efield(field, x, k₀))) × BEAST.n), receivecurrents.functionspace)
    return 0.5 * udot(incident_field, receivecurrents.excitations)
end
function transmission(field, receivecurrents::MagneticSurfaceCurrentDensity{R,S}, k₀::Number) where {R<:Real,S}
    incident_field = BEAST.assemble(((BEAST.n × (x -> hfield(field, x, k₀))) × BEAST.n), receivecurrents.functionspace)
    return -0.5 * udot(incident_field, receivecurrents.excitations)
end
function transmission(field, receivecurrents::ElmagSurfaceCurrentDensity{R,S}, k₀::Number) where {R<:Real,S}
    incident_efield = BEAST.assemble(((BEAST.n × (x -> efield(field, x, k₀))) × BEAST.n), receivecurrents.functionspace)
    incident_hfield = BEAST.assemble(((BEAST.n × (x -> hfield(field, x, k₀))) × BEAST.n), receivecurrents.functionspace)
    return 0.5 *
           (udot(incident_efield, receivecurrents.electricexcitations) - udot(incident_hfield, receivecurrents.magneticexcitations))
end


#### Transformations
function rotate(currents::AbstractSurfaceCurrentDensity, χ::Number, θ::Number, ϕ::Number)
    newcurrents = deepcopy(currents)
    if abs(χ) > 1e-16
        CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, 0.0, χ])
    end
    if abs(θ) > 1e-16
        CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, θ, 0.0])
    end
    if abs(ϕ) > 1e-16
        CompScienceMeshes.rotate!(newcurrents.functionspace.geo, [0.0, 0.0, ϕ])
    end
    return newcurrents
end
function rotate(currents::AbstractSurfaceCurrentDensity; χ=0.0, θ=0.0, ϕ=0.0)
    return rotate(currents, χ, θ, ϕ)
end


function shiftrepresentation(currents::AbstractSurfaceCurrentDensity, R)
    newcurrents = deepcopy(currents)
    CompScienceMeshes.translate!(newcurrents.functionspace.geo, -R)
    return newcurrents
end


