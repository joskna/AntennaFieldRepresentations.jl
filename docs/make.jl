using Documenter, AntennaFieldRepresentations

makedocs(;
    modules=[AntennaFieldRepresentations],
    authors="Josef Knapp <josef.knapp@tum.de>, Danijel Jukic, and Simon B. Adrian",
    sitename="AntennaFieldRepresentations.jl",
    remotes= nothing,
    pages=[
        "Home" => "index.md",
        "Antenna Field Representations" => Any[
            "General Field Representation Interface" => "fieldrepresentations.md",
            "Dipole Representations" => "dipoles/dipole_representations.md",
            "Spherical Representations" => "spherical/spherical_representations.md",
            "Plane Wave Representations" => "planewaves/planewave_representations.md",
            "Equivalent Surface Currents" => "surfacecurrents/surface_representations.md",
        ],
        "Field Samplings"  => "fieldsamplings.md", 
        "OperationMaps" => "operationmaps.md",
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
