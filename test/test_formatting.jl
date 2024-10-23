using JuliaFormatter
pkgpath = pkgdir(AntennaFieldRepresentations)   # path of this package including name
@test format(joinpath(pkgpath, "src", "AntennaFieldRepresentations.jl"), overwrite = false)

@test format(
    joinpath(pkgpath, "src", "SphericalVectorModeFields", "modefunctions.jl"),
    overwrite = false,
)
@test format(
    joinpath(pkgpath, "src", "SphericalVectorModeFields", "modeoperations.jl"),
    overwrite = false,
)
@test format(
    joinpath(pkgpath, "src", "SphericalVectorModeFields", "normalizedlegendre.jl"),
    overwrite = false,
)

@test format(joinpath(pkgpath, "src", "DipoleInteractions", "dipole.jl"), overwrite = false)

@test format(
    joinpath(pkgpath, "src", "PlaneWaveRepresentations", "planewaves.jl"),
    overwrite = false,
)
@test format(
    joinpath(pkgpath, "src", "PlaneWaveRepresentations", "interpolation.jl"),
    overwrite = false,
)

# @test format(joinpath(pkgpath,"src", "Conversions", "DipolePlaneWave.jl"), overwrite=false)
# @test format(joinpath(pkgpath,"src", "Conversions", "DipoleSpherical.jl"), overwrite=false)
# @test format(joinpath(pkgpath,"src", "Conversions", "SphericalPlaneWave.jl"), overwrite=false)
# @test format(joinpath(pkgpath,"src", "Conversions", "SurfaceCurrentPlaneWave.jl"), overwrite=false)
# @test format(joinpath(pkgpath,"src", "Conversions", "SurfaceCurrentSpherical.jl"), overwrite=false)

@test format(
    joinpath(pkgpath, "src", "SurfaceCurrentDensities", "currentrepresentations.jl"),
    overwrite = false,
)

@test format(
    joinpath(pkgpath, "src", "FastSphericalVectorModeTransformations", "wacker.jl"),
    overwrite = false,
)

@test format(joinpath(pkgpath, "src", "MLFMMMatrix", "beastglue.jl"), overwrite = false)
@test format(joinpath(pkgpath, "src", "MLFMMMatrix", "MLFMMMatrix.jl"), overwrite = false)
@test format(
    joinpath(pkgpath, "src", "MLFMMTree", "AbstractMLFMMTree.jl"),
    overwrite = false,
)
@test format(joinpath(pkgpath, "src", "MLFMMTree", "MLFMMTree.jl"), overwrite = false)
@test format(joinpath(pkgpath, "src", "MLFMMTree", "MLFMMTrees.jl"), overwrite = false)
@test format(joinpath(pkgpath, "src", "MLFMMTree", "PBMLFMMTree.jl"), overwrite = false)


