using Documenter, AntennaFieldRepresentations

makedocs(;
    modules = [AntennaFieldRepresentations],
    authors = "Josef Knapp <josef.knapp@tum.de>, Danijel Jukic, and Simon B. Adrian",
    sitename = "AntennaFieldRepresentations.jl",
    remotes = nothing,
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Antenna Field Representations" => Any[
            "General Field Representation Interface"=>"fieldrepresentations.md",
            "Dipole Representations"=>"dipoles/dipole_representations.md",
            "Spherical Mode Representations"=>"spherical/spherical_representations.md",
            "Plane Wave Representations"=>"planewaves/planewave_representations.md",
            "Equivalent Surface Currents"=>"surfacecurrents/surface_representations.md",
        ],
        "Field Samplings" => Any[
            "General Field Sampling Interface"=>"sampling/fieldsamplings.md",
            "Irregularly Distributed Field Samplings"=>"sampling/irregularsampling.md",
            "Spherical Field Samplings"=>"sampling/sphericalsampling.md",
        ],
        "Operations" => Any[
            "OperationMaps"=>"operations/operationmaps.md",
            "Coordinate Transformations"=>"operations/coordinate_trafos.md",
            "Conversion into Other Field Representations"=>"operations/conversions.md",
            "Transmission"=>"operations/transmission.md",
            "Interpolation"=>"operations/interpolation.md",
        ],
        "Examples" => "examples.md",
        "Electromagnetic Theory" => Any[
            "Dipoles"=>"dipoles/dipole_theory.md",
            "Spherical Vector Wave Expansion"=>"spherical/spherical_theory.md",
            "Plane Wave Expansion"=>"planewaves/planewave_theory.md",
            "Equivalent Surface Currents"=>"surfacecurrents/surface_theory.md",
            "Fast Algorithms"=>Any[
                "Spherical Wacker Algorithm" => "spherical/spherical_fast.md"
                "Multilevel Fast Multipole Method" => "mlfmm/mlfmm_theory.md"
            ],
        ],
        "API" => "api.md",
    ],
)
deploydocs(repo = "github.com/joskna/AntennaFieldRepresentations.jl.git")
