using JuliaFormatter
pkgpath = pkgdir(AntennaFieldRepresentations)   # path of this package including name
@test format(joinpath(pkgpath, "."), overwrite = false)
