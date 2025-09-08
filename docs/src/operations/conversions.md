# Conversions Between Field Representations
We can convert an `AntennaRepresentation` into an `AntennaRepresentation` of a different type by the [`changerepresentation(Tnew::Type{<:AntennaFieldRepresentation}, aut_field::AntennaFieldRepresentation; kwargs...)`](@ref) command.

Refer to the following table to check if a certain original representation can be converted into a specific target representation by `AntennaFieldRepresentations`:

| Original Representation  |  |  | Target Representation|  | | | | |
| :------------- | :--------: | :------------: |:------------:  | :------------: |  :------------: |  :------------: |  :------------: | :------------: |
|    |            | [`SphericalWaveExpansion{Radiated}`](@ref)| [`SphericalWaveExpansion{Incident}`](@ref) | [`SphericalWaveExpansion{Absorbed}`](@ref) | [`PlaneWaveExpansion{Radiated}`](@ref)| [`PlaneWaveExpansion{Incident}`](@ref) | [`PlaneWaveExpansion{Absorbed}`](@ref) | [`MLFMMSource`](@ref) |
| [`DipoleArray{Radiated}`](@ref)  |            | ✅ | ✅| |✅ |✅ | | ✅ |
| [`DipoleArray{Incident}`](@ref)  |            |  | ✅| | |✅ | |  |
| [`DipoleArray{Absorbed}`](@ref)  |            |  | |✅ | | |✅ |  |
| [`SurfaceCurrentDensity{Radiated}`](@ref)  |            | ✅ | ✅| |✅ |✅ | | ✅ |
| [`SurfaceCurrentDensity{Incident}`](@ref)  |            |  | ✅| | |✅ | |  |
| [`SurfaceCurrentDensity{Absorbed}`](@ref)  |            |  | |✅ | | |✅ |  |
| [`SphericalWaveExpansion{Radiated}`](@ref)  |            | ✅ | 💨[^1]| |✅|💨[^1] | | |
| [`SphericalWaveExpansion{Incident}`](@ref)  |            |  |✅ | | |✅ | |  |
| [`SphericalWaveExpansion{Absorbed}`](@ref)  |            |  | |✅ | | |✅ |  |
| [`SphericalWaveExpansion{Radiated}`](@ref)  |            | ✅ | 💨[^1]| |✅|💨[^1] | | |
| [`SphericalWaveExpansion{Incident}`](@ref)  |            |  |✅ | | |✅ | |  |
| [`SphericalWaveExpansion{Absorbed}`](@ref)  |            |  | |✅ | | |✅ |  |

[^1]: Cannot be converted into the target representation in the same coordinate system. The coordinate system of the target representation must be translated with repsect to the origina coordinate system. The [`transfer`] method accomplishes this task.

Check out the [conversion examples](@ref hertzarrayexample) in the [examples section](@ref examplessection) for details on the [`changerepresentation`](@ref) command


## ChangeRepresentationMaps
