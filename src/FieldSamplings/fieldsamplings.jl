"""
    ProbeAntenna{A::AntennaFieldRepresentation}

Wrapper around an `AntennaFieldRepresentation` to indicate that it is used as a probe.
"""
struct ProbeAntenna{A <: AntennaFieldRepresentation}
    aut_field::A 
    probesize::Real
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

abstract type FieldSampling{C} end
Base.length(fs::FieldSampling) = prod(size(fs))

"""
    IrregularFieldSampling{P<:ProbeAntenna, T<: Real, C<: Complex}

Field sampling with arbitrary probes at irregularly distributed measurement positions 
"""
struct IrregularFieldSampling{P <: ProbeAntenna, T<: Real, C<: Complex} <: FieldSampling{C}
    positions::Array{SVector{3,T}}
    eulerangles::Array{Tuple{T,T,T}}
    probeIDs::Array{<:Integer}
    probes::Vector{P}
    S21values::Array{C}
    
    function IrregularFieldSampling(
        positions::Array{SVector{3,T}}, 
        eulerangles::Array{Tuple{T,T,T}},
        probeIDs::Array{<:Integer},
        probes::Vector{P},
        S21values::Array{C}) where {P, T, C}

        (size(positions) != size(eulerangles)) && throw(DimensionMismatch("Input dimensions of position list and rotation angle list do not match."))
        (size(positions) != size(probeIDs)) && throw(DimensionMismatch("Input dimensions of position list and probeID list do not match."))
        (size(positions) != size(S21values)) && throw(DimensionMismatch("Input dimensions of position list and S₂₁-list do not match."))
        return new{P, T, C}(positions, eulerangles, probeIDs, probes, S21values)        
    end
end
Base.size(fs::IrregularFieldSampling) = length(fs.S21values)
Base.getindex(fs::IrregularFieldSampling, i) = getindex(fs.S21values, i)
Base.setindex!(fs::IrregularFieldSampling, i, v) = setindex!(fs.S21values, i, v)
function Base.similar(fs::IrregularFieldSampling)
    return IrregularFieldSampling(fs.positions, fs.eulerangles, fs.probeIDs, fs.probes, similar(fs.S21values))
end
function asvector(fs::IrregularFieldSampling)
    return vec(fs.S21values)
end
function IrregularFieldSampling(
    positions::Array{SVector{3,T}}, 
    eulerangles::Array{Tuple{T,T,T}},
    probeIDs::Array{<:Integer},
    probes::Vector{P}) where{P<:ProbeAntenna, T<: Real}

    S21values=zeros(Complex{T}, size(positions))
    return IrregularFieldSampling(positions, eulerangles, probeIDs, probes, S21values)
end

"""
    IrregularFieldSampling(positions::Array{V}, eulerangles::Array{D}, probeIDs::Array{<:Integer}, probes::Vector{P}) where{V, D, P<:ProbeAntenna}

Returns an `IrregularFieldSampling` struct for sampling an `AntennaFieldRepresentation` with arbitrary probes at irregularly distributed measurement positions.

# Inputs:
- `positions::Array{V}` : An array of 3D position vectors. Julia must be able to `convert` the type `V` into an `SVector{3,T<:Real}`.
- `eulerangles::Array{D}` : An array of 3-Tuples denoting the Euler angles `ϑ`, `φ`, and `χ` for rotating the probe antenna at each sample position. Julia must be able to `convert` the type `D` into a `Tuple{T,T,T}}`.
- `probeIDs::Array{<:Integer}` : Array of probeIDs (= indices) to assign a probe from the list `probes` to each sample.
- `probes::Vector{P<:ProbeAntenna}`: Vector of probe antennas which occur in the field sampling 
"""
function IrregularFieldSampling(
    positions::Array{V}, 
    eulerangles::Array{D},
    probeIDs::Array{<:Integer},
    probes::Vector{P}) where{V, D, P<:ProbeAntenna}
    T=eltype(V)
    return IrregularFieldSampling(convert.(SVector{3,T}, positions), convert.(Tuple{T,T,T}, Tuple.(eulerangles)), probeIDs, probes)
end

"""
    EfieldSampling(positions::Vector{V})

Returns an `IrregularFieldSampling` struct for sampling the electric field of an `AntennaFieldRepresentation` at given positions with Hertzian dipole probes.

For an input vector `positions` containing `N` position entries, the resulting `IrregularFieldSampling` provides the `S₂₁`-measurement signal in a `N × 3` matrix, 
where the coloumns correspond to the x-, y-, and z-component of the E-field, respectively.

# Input:
- `positions::Vector{V}` : A vector of 3D position vectors. Julia must be able to `convert` the type `V` into an `SVector{3,T<:Real}`. 
"""
function EfieldSampling(positions::Vector{V}) where{V}
    T= eltype(V)
    pos1= convert.(SVector{3, T}, positions)
    pos= [pos1 pos1 pos1]

    tuple1= convert(Tuple{T,T,T}, Tuple([0,0,0])) # for z-oriented dipole
    tuple2= convert(Tuple{T,T,T}, Tuple([pi/2,0,0])) # for x-oriented dipole
    tuple3= convert(Tuple{T,T,T}, Tuple([pi/2,pi/2,0])) # for y-oriented dipole
    eulerangles= Array{Tuple{T,T,T}}(undef, length(positions), 3)
    for k in size(eulerangles, 1) 
        eulerangles[k,1] =  tuple2
        eulerangles[k,2] =  tuple3
        eulerangles[k,3] =  tuple1
    end
    
    probeIDs= ones(Int, size(pos))
    probes=[ProbeAntenna(HertzianArray([T.([0, 0, 0])], [Complex{T}.([0, 0, 1])], [Complex{T}.(2)], Complex{T}(0.0)), 0.0)]
    
    return IrregularFieldSampling(pos, eulerangles, probeIDs, probes)
end

"""
    HfieldSampling(positions::Vector{V})

Returns an `IrregularFieldSampling` struct for sampling the magnetic field of an `AntennaFieldRepresentation` at given positions with Fitzgerald dipole probes.

For an input vector `positions` containing `N` position entries, the resulting `IrregularFieldSampling` provides the `S₂₁`-measurement signal in a `N × 3` matrix, 
where the coloumns correspond to the x-, y-, and z-component of the H-field, respectively.

# Input:
- `positions::Vector{V}` : A vector of 3D position vectors. Julia must be able to `convert` the type `V` into an `SVector{3,T<:Real}`. 
"""
function HfieldSampling(positions::Vector{V}) where{V}
    T= eltype(V)
    pos1= convert.(SVector{3, T}, positions)
    pos= [pos1 pos1 pos1]

    tuple1= convert(Tuple{T,T,T}, Tuple([0,0,0])) # for z-oriented dipole
    tuple2= convert(Tuple{T,T,T}, Tuple([pi/2,0,0])) # for x-oriented dipole
    tuple3= convert(Tuple{T,T,T}, Tuple([pi/2,pi/2,0])) # for y-oriented dipole
    eulerangles= Array{Tuple{T,T,T}}(undef, length(positions), 3)
    for k in size(eulerangles, 1) 
        eulerangles[k,1] =  tuple2
        eulerangles[k,2] =  tuple3
        eulerangles[k,3] =  tuple1
    end
    
    probeIDs= ones(Int, size(pos))
    probes=[ProbeAntenna(FitzgeraldArray([T.([0, 0, 0])], [Complex{T}.([0, 0, 1])], [Complex{T}.(2)], Complex{T}(0.0)), 0.0)]
    
    return IrregularFieldSampling(pos, eulerangles, probeIDs, probes)
end

"""
    RegularSphericalFieldSampling{Y<:SphereSamplingStrategy, S<:AbstractSphericalCoefficients, C<:Complex}
    
Field sampling on spherical measurement surface with measurement positions distributed according to a `SphereSamplingStrategy`
"""
struct RegularSphericalFieldSampling{Y<:SphereSamplingStrategy,S<:AbstractSphericalCoefficients,C<:Complex} <: FieldSampling{C}
    samplingstrategy::Y
    incidentcoefficients::S
    S21values::Matrix{C}
end
Base.size(fs::RegularSphericalFieldSampling) = length(fs.S21values)
Base.getindex(fs::RegularSphericalFieldSampling, i) = getindex(fs.S21values, i)
Base.setindex!(fs::RegularSphericalFieldSampling, i, v) = setindex!(fs.S21values, i, v)
function Base.similar(fs::RegularSphericalFieldSampling)
    return RegularSphericalFieldSampling(fs.samplingstrategy, fs.incidentcoefficients, similar(fs.S21values))
end
function asvector(fs::RegularSphericalFieldSampling)
    return vec(fs.S21values)
end

"""
    FarFieldSampling{Y<:SphereSamplingStrategy, C<:Complex}

Field sampling of the far-field ϑ- and ϕ-components with measurement positions distributed according to a `SphereSamplingStrategy`
"""
struct FarFieldSampling{Y<:SphereSamplingStrategy, C<:Complex} <: FieldSampling{C}
    samplingstrategy::Y
    Eθϕ::Matrix{C}
end
Base.size(fs::FarFieldSampling) = length(fs.Eθϕ)
Base.getindex(fs::FarFieldSampling, i) = getindex(fs.Eθϕ, i)
Base.setindex!(fs::FarFieldSampling, i, v) = setindex!(fs.Eθϕ, i, v)
function Base.similar(fs::FarFieldSampling)
    return FarFieldSampling(fs.samplingstrategy, similar(fs.Eθϕ))
end
function asvector(fs::FarFieldSampling)
    return vec(fs.Eθϕ)
end

"""
    FarFieldSampling{Y<:SphericalSamplingStrategy, C<:Complex}

Field sampling of the far-field ϑ- and ϕ-components with measurement positions irregularly distributed on the far-field sphere
"""
struct IrregularFarFieldSampling{C<:Complex, T1<:Number, T2<:Real} <: FieldSampling{C}
    θϕ::Vector{Tuple{T1,T2}}
    Eθϕ::Vector{C}
end
Base.size(fs::IrregularFarFieldSampling) = length(fs.Eθϕ)
Base.getindex(fs::IrregularFarFieldSampling, i) = getindex(fs.Eθϕ, i)
Base.setindex!(fs::IrregularFarFieldSampling, i, v) = setindex!(fs.Eθϕ, i, v)
function Base.similar(fs::IrregularFarFieldSampling)
    return IrregularFarFieldSampling(fs.θϕ, similar(fs.Eθϕ))
end
function asvector(fs::IrregularFarFieldSampling)
    return vec([fs.Eθ; fs.Eϕ])
end
function _eθ(fs::IrregularFarFieldSampling)
    return view(fs.Eθϕ, :, 1:(length(fs) ÷ 2))
 end
 function _eϕ(fs::IrregularFarFieldSampling)
    return view(fs.Eθϕ, :, (length(fs) ÷ 2 + 1) : length(fs))
 end

 """
    SphericalFieldSampling

Field sampling on regularly distributed measurement positions on a spherical surface
"""
struct SphericalFieldSampling{Y<:SphereSamplingStrategy, A<:AbstractSphericalCoefficients, C<:Complex} <: FieldSampling{C}
    incidentcoefficients::A
    samplingstrategy::Y
    S21values::Array{C}   
end
function SphericalFieldSampling(samplingstrategy::Y, incidentcoefficiens::A) where{C, Y <: SphereSamplingStrategy, A <:AbstractSphericalCoefficients{C}}
    samplecountθ, samplecountϕ=_countsamples(samplingstrategy)
    S21values= zeros(C, samplecountθ , samplecountϕ, 2)
    return SphericalFieldSampling{Y, A, C}(incidentcoefficiens, samplingstrategy, S21values)
end
function asvector(x::Array{C}) where{C}
    return vec(x)
end
Base.size(fs::SphericalFieldSampling) = (length(fs.S21values), )
Base.getindex(fs::SphericalFieldSampling, i) = getindex(fs.S21values, i)
Base.setindex!(fs::SphericalFieldSampling, i, v) = setindex!(fs.S21values, i, v)
function Base.similar(fs::SphericalFieldSampling)
    return SphericalFieldSampling(deepcopy(fs.incidentcoefficients), fs.samplingstrategy, similar(fs.S21values))
end