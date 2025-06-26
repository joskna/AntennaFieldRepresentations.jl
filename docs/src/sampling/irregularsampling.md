# [Irregularly Distributed Field Samplings](@id irregularsampling)
The most general [`FieldSampling`](@ref fieldsampling) in `AntennaFieldRepresentations.jl` is the `IrregularFieldSampling`. 
It allows to sample the antenna field at arbitrary positions with arbitrary probe antennas.
The user must specify
- a list of measurement positions 
- a vector of [probe antennas](@ref probeantenna) which are used in the measurement scenario
- a list of probe IDs to assign the utilized probe to each measurement position
- a list of Euler angles ``\vartheta``, ``\varphi``, ``\chi`` to define the probe rotation at each measurement position

## Constructors of `IrregularFieldSampling`
To create an [`IrregularFieldSampling`](@ref), use the following constructor:
```julia
function IrregularFieldSampling(
    positions,
    eulerangles,
    probeIDs,
    probes
)
```

with the following inputs:
- `positions::Array{V}` : An array of 3D position vectors. Julia must be able to `convert` the type `V` into an `SVector{3,T<:Real}`.
- `eulerangles::Array{D}` : An array of 3-Tuples denoting the Euler angles `ϑ`, `φ`, and `χ` for rotating the probe antenna at each sample position. Julia must be able to `convert` the type `D` into a `Tuple{T,T,T}}`.
- `probeIDs::Array{<:Integer}` : Array of probeIDs (= positive indices) to assign a probe from the list `probes` to each sample.
- `probes::Vector{P<:ProbeAntenna}`: Vector of probe antennas which occur in the field sampling 

The arrays `positions`, `eulerangles`, and `probIDs` must be of the same size (i.e. same number of dimenions and same number of entries along the respective dimensions). The length of the vector `probes` must not be smaller than the largest element of `probeIDs`.
This way, every measurement sample is assigned to a position in space, a probe antenna from the probe antennas listed in `probes` and an orientation of the probe antenna.

## E-Field and H-Field Requests
Since unfiltered sampling of the electric or magnetic field of an [`AntennaFieldRepresentation`](@ref fieldrepresentation) at given measurement locations is a common use case, a convenience constructor is provided by `AntennaFieldRepresentations.jl` for these cases. The methods `EfieldSampling(positions::Vector{V})` and `HfieldSampling(positions::Vector{V})` return `IrregularFieldSampling`s which correspond to sampling the field of an `AntennaFieldRepresentation` at given positions with [Hertzian dipole probes](@ref receivedsignaldipoles) or [Fitzgerald dipole probes](@ref receivedsignaldipoles) respectively.