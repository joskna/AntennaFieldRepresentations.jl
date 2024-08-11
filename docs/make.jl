using Documenter, AntennaFieldRepresentations

makedocs(;
    modules=[AntennaFieldRepresentations],
    authors="Josef Knapp <josef.knapp@tum.de>, Danijel Jukic, and Simon B. Adrian",
    sitename="AntennaFieldRepresentations.jl",
    remotes= nothing,
    pages=[
        "Home" => "index.md",
        "Antenna Field Representations"  => "fieldrepresentations.md", 
        "Operations" => Any[
            "Field Evaluations" => "operations/fields.md",
            "Coordinate Transformations" => "operations/coordinate_trafos.md",
            "Conversions" => "operations/conversions.md",
            "Transmissions" => "operations/interactions.md"
        ],
        "Examples" => "examples.md",
        "Electromagnetic Theory"=> Any[
            "Dipoles" => "dipoles/dipole_theory.md",
            "Spherical Vector Wave Expansion" => "spherical/spherical_theory.md",
            "Plane Wave Expansion" => "planewaves/planewave_theory.md",
            "Equivalent Surface Currents" => "surfacecurrents/surface_theory.md",
            "Fast Algorithms" => Any[
                "Spherical Wacker Algorithm" => "spherical/spherical_fast.md"
                "Multilevel Fast Multipole Method" => "mlfmm/mlfmm_theory.md"
            ]
        ], 
        "API" => "api.md",
        
    ],
)
deploydocs(
    repo="github.com/joskna/AntennaFieldRepresentations.jl.git"
)
