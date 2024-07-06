module AntennaFieldRepresentations
using StaticArrays
using LinearAlgebra: norm, BLAS.dotu, ⋅, cross, rank, mul!

using SpecialFunctions: hankelh1, hankelh2, besselj

using WignerD: wignerd

using FastGaussQuadrature: gausslegendre
# using FourierTools: resample

using BEAST
using CompScienceMeshes
using ClusterTrees

using FFTW
# using Memoization

export transmission, rotate, translate, shiftrepresentation
export efield,
    hfield,
    ehfield,
    farfield,
    converttype,
    convertrepresentation,
    elementtype,
    reciprocaltype
export udot, sqrtdot, cdist, rot_mat_zyz
export c₀, ε₀, μ₀, Z₀
export AntennaFieldRepresentation

# using PlotlyJS
# import PlotlyJS: plot

using LinearMaps
# using IterativeSolvers

const c₀ = 299792458 # m/s
const ε₀ = 8.854187812813e-12 # F/m
const μ₀ = 1.2566370621219e-6 # N/A²
const Z₀ = sqrt(μ₀ / ε₀) # Ω


abstract type AntennaFieldRepresentation end


"""
    transmission(sender<:FieldRepresentation, receiver<:FieldRepresentation, k₀::Number) -> val::Complex

Return interaction between two field representations ( or source representations)
"""
function transmission(sender, receiver::T, k₀::Number) where {T<:AntennaFieldRepresentation}
    return transmission(convertrepresentation(reciprocaltype(T), sender, k₀), receiver, k₀)
end
# function transmission(sender::Array{<:AntennaFieldRepresentation}, receiver::T, k₀::Number) where{T<:AntennaFieldRepresentation}
#     return transmission(convertrepresentation(reciprocaltype(T), sender), receiver, k₀)
# end


"""
    elementtype(fieldobject<:FieldRepresentation)

Return base type of elements representing the fields in `fieldobject`. Default is ComplexF64
"""
function elementtype end

"""
    rotate(fieldobject<:FieldRepresentation, χ::Number, θ::Number, ϕ::Number) -> rotatedfieldobject<:Fieldrepresentation
    rotate(fieldobject<:FieldRepresentation; [χ=0.0],  [θ=0.0], [ϕ=0.0]) -> rotatedfieldobject<:Fieldrepresentation

Rotate a field representation around the Euler angles `χ`, `θ`, `ϕ`.
"""
function rotate(fieldrepresentation; χ = 0.0, θ = 0.0, ϕ = 0.0)
    return rotate(fieldrepresentation, χ, θ, ϕ)
end

"""
    translate(fieldobject<:FieldRepresentation, R, k₀::Number) -> newobject<:Fieldrepresentation

Translate radiating field representation into incident field representation in translated coordinate system with new origin at R. 
"""
function translate end

"""
    shiftrepresentation(fieldobject<:FieldRepresentation, R, k₀::Number) -> newobject<:Fieldrepresentation

Translate field representation into field representation of same type in translated coordinate system with new origin at R. 
"""
function shiftrepresentation end


"""
    reciprocaltype(type::Type{<:AntennaFieldRepresentation}) -> Type{<:AntennaFieldRepresentation}
Return type of the field representation which is reciprocal to the input type.    
"""
function reciprocaltype(type::Type{<:AntennaFieldRepresentation}) end


"""
    efield(fieldobject<:FieldRepresentation, R, k₀::Number) -> [Ex;Ey;Ez]
    efield(Vector{<:FieldRepresentation}, R, k₀::Number) -> [Ex;Ey;Ez]

Return the E-field vector (in cartesian coordinates) of the field representation at location R.
"""
function efield(fieldrepresentation, Rvec::AbstractVector{<:Number}, k₀::Number)
    return efield(fieldrepresentation, [Rvec], k₀)[1]
end
function efield(
    fields::Array{<:AntennaFieldRepresentation},
    R::AbstractVector{<:Number},
    k₀::Real,
)
    return sum(efield(field, R, k₀) for field in fields)
end
function efield(fieldrepresentation, Rvecs::Array{<:AbstractVector{<:Number}}, k₀::Number)
    return [efield(fieldrepresentation, Rvec, k₀) for Rvec in Rvecs]
end


"""
    hfield(fieldobject<:FieldRepresentation, R, k₀::Number) -> [Hx;Hy;Hz]
    hfield(Vector{<:FieldRepresentation}, R, k₀::Number) -> [Hx;Hy;Hz]

Return the H-field vector (in cartesian coordinates) of the field representation at location R.
"""
function hfield(fieldrepresentation, Rvec::AbstractVector{<:Number}, k₀::Number)
    return hfield(fieldrepresentation, [Rvec], k₀)[1]
end
function hfield(
    fields::Array{<:AntennaFieldRepresentation},
    R::AbstractVector{<:Number},
    k₀::Real,
)
    return sum(hfield(field, R, k₀) for field in fields)
end
function hfield(fieldrepresentation, Rvecs::Array{<:AbstractVector{<:Number}}, k₀::Number)
    return [hfield(fieldrepresentation, Rvec, k₀) for Rvec in Rvecs]
end

"""
    ehfield(fieldobject<:FieldRepresentation, R, k₀::Number) -> [Ex;Ey;Ez], [Hx;Hy;Hz]
    ehfield(Vector{<:FieldRepresentation}, R, k₀::Number) -> [Ex;Ey;Ez], [Hx;Hy;Hz]

Return the E-field and H-field vector (in cartesian coordinates) of the field representation at location R.
"""
function ehfield(fieldrepresentation, Rvec::AbstractVector{<:Number}, k₀::Number)
    return ehfield(fieldrepresentation, [Rvec], k₀)[1]
end
function ehfield(
    fields::Array{<:AntennaFieldRepresentation},
    R::AbstractVector{<:Number},
    k₀::Real,
)
    return sum.(ehfield(field, R, k₀) for field in fields)
end
# function ehfield(fieldrepresentation, Rvecs::Array{<:AbstractVector{<:Number}}, k₀::Number)
#     return [ehfield(fieldrepresentation, Rvec, k₀) for Rvec in Rvecs]
# end

"""
    farfield(fieldrepresentation, θ, ϕ, k₀::Number) -> Tuple(Complex64,2)   

Return Eθ,Eϕ-far-field tuple for  of radiating field representation into direction(s) specified by θ and ϕ.
"""
# function farfield(fieldrepresentation, θ::Number, ϕ::Number, k₀::Number)
#     Eθ, Eϕ = farfield(fieldrepresentation, [θ], [ϕ], k₀)
#     return Eθ[1, 1], Eϕ[1, 1]
# end
function farfield(
    fieldrepresentation,
    θvec::Vector{<:Number},
    ϕvec::Vector{<:Number},
    k₀::Number,
)
    Etuples = [farfield(fieldrepresentation, θ, ϕ, k₀) for θ in θvec, ϕ in ϕvec]
    return [e[1] for e in Etuples], [e[2] for e in Etuples]
end
function farfield(fieldrepresentation, Jθ::Integer, Jϕ::Integer, k₀::Number)
    ϕvec = collect(range(0; stop = 2π, length = Jϕ + 1))[1:(end-1)]
    θvec = collect(range(0; stop = 2π, length = Jθ + 1))[1:Integer(floor((Jθ + 1) / 2))]
    return farfield(fieldrepresentation, θvec, ϕvec, k₀)
end
function farfield(
    fieldrepresentation,
    θϕtuples::Vector{Tuple{<:Number,<:Number}},
    k₀::Number,
)
    Etuples =
        [farfield(fieldrepresentation, θϕtuple[1], θϕtuple[2], k₀) for θϕtuple in θϕtuples]
    return [e[1] for e in Etuples], [e[2] for e in Etuples]
end


"""
    converttype(newType<:Type{FieldRepresentation}, fieldobject<: FieldRepresentation)

Convert field representation into other type with same datastructure.
Field avaluations of the respective old and new field representations will differ.
"""
function converttype end


"""
    convertrepresentation(newType<:Type{FieldRepresentation}, fieldobject<: FieldRepresentation, args)

 Convert field representation into other type with different datastructure.
 
 Field evaluations of the respective old and new field representations should be approximately equal in a common region of convergence.
"""
function convertrepresentation end

"""
    dotu(u::Array{<:Number,1}, v::Array{<:Number,1}) -> Number

Return u⋅v without complex conjugation of u
"""
function udot(u::Array{<:Real,1}, v::Array{<:Real,1})
    return u ⋅ v
end
function udot(u::Array{<:Complex,1}, v::Array{<:Complex,1})
    return dotu(length(u), u, 1, v, 1)
end
function udot(u::Array{<:Complex,1}, v::Array{<:Real,1})
    return dotu(length(u), u, 1, complex(v), 1)
end
function udot(u::Array{<:Real,1}, v::Array{<:Complex,1})
    return dotu(length(u), complex(u), 1, v, 1)
end
function udot(u, v)
    # return sum(u .* v)
    sumval = zero(promote_type(eltype(u), eltype(v)))
    for k in eachindex(u)
        sumval += u[k] .* v[k]
    end
    return sumval
end



"""
    sqrtdot(u::Array{<:Number,1}, v::Array{<:Number,1}) -> Number

Return √(u⋅v) without complex conjugation of u
"""
function sqrtdot(u, v)
    return sqrt(udot(u, v))
end


"""
    cdist(R::Array{<:Number,1}) -> Number

Return √(R⋅R) without complex conjugation of R in the dot product.
Corresponds to the norm operator for real valued vectors
"""
function cdist(R::Array{<:Real,1})
    return norm(R)
end
function cdist(R::StaticArraysCore.SVector{3,<:Real})
    return norm(R)
end
function cdist(R::StaticArraysCore.SVector{3,<:Complex})
    return sqrtdot(R, R)
end
function cdist(R::Array{<:Number,1})
    return sqrtdot(R, R)
end
function cdist(R)
    return sqrtdot(R, R)
end


"""
    rot_mat_zyz(χ::Number, θ::Number, ϕ::Number)-> R::Matrix{eltype(χ, θ, ϕ)}

Return rotation matrix for cartesian vector components for a rotation around Euler angles `χ`, `θ`, `ϕ`.
"""
function rot_mat_zyz(χ::Number, θ::Number, ϕ::Number)
    s1, c1 = sincos(χ)
    s2, c2 = sincos(θ)
    s3, c3 = sincos(ϕ)


    return [
        c1*c2*c3-s1*s3 -s1*c2*c3-c1*s3 s2*c3
        c1*c2*s3+s1*c3 -s1*c2*s3+c1*c3 s2*s3
        -c1*s2 s1*s2 c2
    ]
end

function rot_mat_zyx(yaw::Number, pitch::Number, roll::Number)

    # compute cos and sin
    sin_yaw, cos_yaw = sincos(yaw)
    sin_pitch, cos_pitch = sincos(pitch)
    sin_roll, cos_roll = sincos(roll)


    # rotation matrix
    R = zeros(Float64, 3, 3)
    R[1, :] = [cos_pitch * cos_yaw, cos_pitch * sin_yaw, -sin_pitch]
    R[2, :] = [
        sin_roll * sin_pitch * cos_yaw - cos_roll * sin_yaw,
        sin_roll * sin_pitch * sin_yaw + cos_roll * cos_yaw,
        sin_roll * cos_pitch,
    ]
    R[3, :] = [
        cos_roll * sin_pitch * cos_yaw + sin_roll * sin_yaw,
        cos_roll * sin_pitch * sin_yaw - sin_roll * cos_yaw,
        cos_roll * cos_pitch,
    ]

    return R

end


## SphericalVectorModeFields
include(joinpath("SphericalVectorModeFields", "normalizedlegendre.jl"))
include(joinpath("SphericalVectorModeFields", "modefunctions.jl"))
include(joinpath("SphericalVectorModeFields", "radialfunctions.jl"))
include(joinpath("SphericalVectorModeFields", "modeoperations.jl"))
export RadiatingSphericalExpansion,
    AbsorbedSphericalExpansion, IncidentSphericalExpansion, AbstractSphericalExpansion
export sℓm_to_j, j_to_sℓm, αtoβ, βtoα

## DipoleInteractions
include("DipoleInteractions/dipole.jl")
export HertzDipole, FitzgeraldDipole, AbstractDipole

## PlaneWaveRepresentations
include(joinpath("PlaneWaveRepresentations", "planewaves.jl"))
include(joinpath("PlaneWaveRepresentations", "interpolation.jl"))
include(joinpath("PlaneWaveRepresentations", "transfer.jl"))
export FarfieldPattern, PlaneWaveSpectrum, PlaneWaveRepresentation, PlaneWave
export samplingrule, resample, revertdirection

## SurfaceCurrentDensities
include(joinpath("SurfaceCurrentDensities", "currentrepresentations.jl"))
export ElectricSurfaceCurrentDensity,
    MagneticSurfaceCurrentDensity, ElmagSurfaceCurrentDensity, AbstractSurfaceCurrentDensity

## FastSphericalVectorModeTransformations
include(joinpath("FastSphericalVectorModeTransformations", "wacker.jl"))
export Wackerforward, Wackerforward_ad, Wacker

## Conversions
include(joinpath("Conversions", "DipoleSpherical.jl"))
include(joinpath("Conversions", "DipolePlaneWave.jl"))
include(joinpath("Conversions", "SphericalPlaneWave.jl"))
include(joinpath("Conversions", "SurfaceCurrentPlaneWave.jl"))
include(joinpath("Conversions", "SurfaceCurrentSpherical.jl"))

## MLFMMTree
include(joinpath("MLFMMTree", "MLFMMTrees.jl"))

#####
include(joinpath("AntennaRepresentations", "antennarepresentations.jl"))
include(joinpath("FieldSamplings", "fieldsamplings.jl"))


##MLFMMMatrix
include(joinpath("MLFMMMatrix", "MLFMMMatrix.jl"))
include(joinpath("MLFMMMatrix", "beastglue.jl"))


## Plotting
# include(joinpath("Plotting", "plotting.jl"))

## SourceReconstruction
include(joinpath("InteractionMatrices", "interactionmatrices.jl"))
export DipoleInteractionMatrix, MLFMMInteractionMatrix
# using ClusterTrees: ClusterTrees
# using StaticArrays: StaticArrays, SVector
# using LinearAlgebra: LinearAlgebra
# using Exceptions: Exceptions



# include(joinpath("MLFMMTree", "AbstractMLFMMTree.jl"))
# include(joinpath("MLFMMTree", "MLFMMTree.jl"))
# export getboundingbox



end
