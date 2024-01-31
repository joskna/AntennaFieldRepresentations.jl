using Documenter, AntennaFieldRepresentations

makedocs(;
    modules=[AntennaFieldRepresentations],
    authors="Josef Knapp <josef.knapp@tum.de>, Danijel Jukic, and Simon B.Adrian",
    sitename="AntennaFieldRepresentations.jl",
    remotes= nothing,
    pages=[
        "Home" => "index.md",
        "Representations" =>Any[
            "Dipoles" => [
                "Theory" =>"dipoles/dipole_theory.md",
                "Examples" => "dipoles/dipole_examples.md"],
            "Spherical Vector Wave Expansion" => Any[
                "Theory" => "spherical/spherical_theory.md",
                "Fast Algorithms" => "spherical/spherical_fast.md",
                "Examples" => "spherical/spherical_examples.md"
            ],
            "Plane Wave Expansion" => 
            Any[
                "Theory" => "planewaves/planewave_theory.md",
                "Examples" => "planewaves/planewave_examples.md"
            ],
            "Equivalent Surface Currents" => 
            Any[
                "Theory" => "surfacecurrents/surface_theory.md",
                "Examples" => "surfacecurrents/surface_examples.md"
            ],
        ],
        "Operations" => Any[
            "Field Evaluations" => "operations/fields.md",
            "Coordinate Transformations" => "operations/coordinate_trafos.md",
            "Conversions" => "operations/conversions.md",
            "Transmissions" => "operations/interactions.md",
        ],
        "Examples" => "examples.md",
        "API" => "api.md",
        
    ],
)
deploydocs(
    repo="github.com/joskna/AntennaFieldRepresentations.jl.git"
)
