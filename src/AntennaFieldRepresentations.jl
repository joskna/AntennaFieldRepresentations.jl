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

# using PlotlyJS
# import PlotlyJS: plot

using LinearMaps
# using IterativeSolvers

const c₀ = 299792458 # m/s
const ε₀ = 8.854187812813e-12 # F/m
const μ₀ = 1.2566370621219e-6 # N/A²
const Z₀ = sqrt(μ₀ / ε₀) # Ω


"""
    udot(u::Array{<:Number,1}, v::Array{<:Number,1}) -> Number

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


include(joinpath("AntennaRepresentations", "antennarepresentations.jl"))
include(joinpath("DipoleInteractions", "dipole.jl"))
include(joinpath("PlaneWaveRepresentations", "planewave.jl"))
include(joinpath("SphericalVectorModeFields", "spherical.jl"))
include(joinpath("SurfaceCurrentDensities", "currentrepresentations.jl"))
include(joinpath("FieldSamplings", "fieldsamplings.jl"))
include(joinpath("OperationMaps", "operationmaps.jl"))
include(joinpath("OperationMaps", "changerepresentationmaps.jl"))

export PropagationType, Radiated, Absorbed, Incident
export AntennaFieldRepresentation
export ElmagType, Electric, Magnetic
export asvector, efield, efield!, hfield, hfield!, ehfield, ehfield!, farfield
export getwavenumber, setwavenumber!, rotate, rotate!, spatialshift, spatialshift!
export DipoleArray, HertzArray, FitzgeraldArray
export SphericalWaveExpansion,
    AbstractSphericalCoefficients, SphericalCoefficients, FirstOrderSphericalCoefficients
export SphereSamplingStrategy,
    PlaneWaveExpansion, RegularθRegularϕSampling, GaussLegendreθRegularϕSampling
export SurfaceCurrentDensity
export FieldSampling, IrregularFieldSampling, ProbeAntenna, EfieldSampling, HfieldSampling
export RegularSphericalFieldSampling
export j_to_sℓm, sℓm_to_j
export getwavenumber, equivalentorder
export changerepresentation, transmit
export ChangeRepresentationMap, TransmitMap
export SphericalFieldSampling, SphericalTransmitMap
export inverse


end
