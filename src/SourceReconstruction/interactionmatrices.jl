
"""
    MLFMMInteractionMatrix{C} <: LinearMaps.LinearMap{C}

Represent interactions between array of dipole sources and array of dipole probes as abstract matrix
# Fields 
- sourcestruct::AntennaFieldRepresentations.MLFMMSource
- receivestruct::AntennaFieldRepresentations.MLFMMReceive
"""
struct MLFMMInteractionMatrix{C} <: LinearMaps.LinearMap{C}
    sourcestruct::AntennaFieldRepresentations.MLFMMSource
    receivestruct::AntennaFieldRepresentations.MLFMMReceive
end
function LinearMaps._unsafe_mul!(y, A::MLFMMInteractionMatrix{C}, x::V) where{C, V<:AbstractVector}
    AntennaFieldRepresentations._forward!(A.sourcestruct, A.receivestruct,x)
    y .= copy(A.receivestruct.bvector)
end
# function LinearMaps._unsafe_mul!(y, B::LinearMaps.AdjointMap{C, MLFMMInteractionMatrix{C}}, x::V) where{C, V<:AbstractVector}
#     A=B.lmap
# A.receivestruct.bvector .= x
#     AntennaFieldRepresentations._adjoint_forward!(A.sourcestruct, A.receivestruct)
#     y .= copy(A.sourcestruct.xvector)
# end
function LinearMaps._unsafe_mul!(y, B::LinearMaps.AdjointMap{C, MLFMMInteractionMatrix{C}}, x::V) where{C, V<:AbstractVector}
    A=B.lmap
    A.receivestruct.bvector .= x
    conj!(A.receivestruct.bvector)
    AntennaFieldRepresentations._transpose_forward!(A.sourcestruct, A.receivestruct)
    conj!(A.sourcestruct.xvector)
    y .= A.sourcestruct.xvector
end
function LinearMaps._unsafe_mul!(y, B::LinearMaps.TransposeMap{C, MLFMMInteractionMatrix{C}}, x::V) where{C, V<:AbstractVector}
    A=B.lmap
    AntennaFieldRepresentations._transpose_forward!(A.sourcestruct, A.receivestruct,x)
    y .= A.sourcestruct.xvector
end

Base.size(A::MLFMMInteractionMatrix) =(length(A.receivestruct.basisfunctionpatterns), length(A.sourcestruct.basisfunctionfarfields))

"""
    DipoleInteractionMatrix <: AbstractArray{ComplexF64,2}

Represent interactions between array of dipole sources and array of dipole probes as abstract matrix
# Fields 
- `sourcedipoles::Array{<:AbstractDipole,1}`
- `probedipoles::Array{<:AbstractDipole,1}` 
- `k₀::Float64`: Wavenumber
"""
# struct DipoleInteractionMatrix{C} <: LinearMaps.LinearMap{C}
struct DipoleInteractionMatrix <: AbstractArray{ComplexF64,2}
    sourcedipoles::Array{<:AbstractDipole,1}
    probedipoles::Array{<:AbstractDipole,1}
    k₀::Float64
end
Base.size(D::DipoleInteractionMatrix) = (length(D.probedipoles), length(D.sourcedipoles))
Base.getindex(D::DipoleInteractionMatrix, i::Int, j::Int) = transmission(D.sourcedipoles[j], D.probedipoles[i], D.k₀)
Base.getindex(D::DipoleInteractionMatrix, i::Number, j::Number) = D[convert(Int, i), convert(Int, j)]
# function LinearMaps._unsafe_mul!(y, A::DipoleInteractionMatrix, x::AbstractVector)
#     Threads.@threads for i in 1:size(A,1)
#         y[i]=sum([transmission(A.sourcedipoles[k], A.probedipoles[i], A.k₀) .* x[k] for k in 1:size(A,2)])
#     end
#     return y
# end
# function LinearMaps._unsafe_mul!(y, B::LinearMaps.AdjointMap{C,DipoleInteractionMatrix{C}}, x::AbstractVector) where{C}
#     A=B.lmap
#     Threads.@threads for i in 1:size(A,1)
#         y[i]=sum([transmission(A.sourcedipoles[i], A.probedipoles[k], A.k₀) .* x[k] for k in 1:size(A,1)])
#     end
#     return y
# end
# function DipoleInteractionMatrix(sourcedipoles::Vector{<:AbstractDipole}, probedipoles::Vector{<:AbstractDipole}, k0::T) where{T<:Real}
#     return DipoleInteractionMatrix{Complex{T}}(sourcedipoles, probedipoles, k0)
# end


struct SurfaceCurrentDipoleInteractionMatrix <: AbstractArray{ComplexF64,2}
    surfacecurrents::AbstractSurfaceCurrentDensity
    probedipoles::Array{<:AbstractDipole,1}
    k₀::Float64
end
Base.size(D::SurfaceCurrentDipoleInteractionMatrix) = (length(D.probedipoles), length(D.surfacecurrents.excitations))
function Base.getindex(D::SurfaceCurrentDipoleInteractionMatrix, i::Int, j::Int)  
    xin=zeros(ComplexF64, size(D,2))
    xin[j]=1
    eqcurrents=typeof(D.surfacecurrents)(D.surfacecurrents.functionspace, xin)
    b=transmission(D.probedipoles[i], eqcurrents, D.k₀)
    return b
end
function Base.getindex(D::SurfaceCurrentDipoleInteractionMatrix, i::UnitRange, j::Int)  
    xin=zeros(ComplexF64, size(D,2))
    xin[j]=1
    eqcurrents=typeof(D.surfacecurrents)(D.surfacecurrents.functionspace, xin)
    b=Vector(ComplexF64,length(i))
    for (index, ii) in enumerate i
        b[index]=transmission(D.probedipoles[ii], eqcurrents, D.k₀)
    end
    
    return b
end
function Base.getindex(D::SurfaceCurrentDipoleInteractionMatrix, i::Int, j::UnitRange)  
    incident_field = 0.5 * BEAST.assemble(((BEAST.n × (x -> efield(D.probedipoles[i], x, D.k₀))) × BEAST.n), D.surfacecurrents.functionspace)
    return incident_field[j]
end
