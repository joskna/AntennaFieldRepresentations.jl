struct MLFMMSource{R<:ResampleMap, Y<: SphereSamplingStrategy, R<:Real}
    tree::MLFMMTree
    expectedaccuracy::R
    wavenumber::R
    nodefarfields::Vector{PlaneWaveExpansion{Radiated,Y,Complex{R}}}
    basisfunctionfarfields::Vector{PlaneWaveExpansion{Radiated,Y,Complex{R}}}
    levelcutoffparameters::Vector{Int}
    levelinterpolators::Vector{I}
    phaseshifttoparent::Array{Matrix{Complex{R}},2}
    nodeisfresh::Vector{Bool}
    leafnodeindices::Vector{Int}
    xvector::Vector{Complex{R}}
    verbose::Bool
    tmpmatrix::Matrix{Complex{R}}
end

struct MLFMMReceive{R<:ResampleMap, Y<: SphereSamplingStrategy, R<:Real}
    tree::MLFMMtree
    # expectedaccuracy::R
    # wavenumber::R
    nodespectra::Vector{PlaneWaveExpansion{Incident,Y,Complex{R}}}
    basisfunctionpatterns::Vector{PlaneWaveExpansion{Radiated,Y,Complex{R}}}
    levelcutoffparameters::Vector{Int}
    levelinterpolators::Vector{A}
    phaseshifttoparent::Array{Matrix{Complex{R}},2}
    transferlist::Vector{Vector{Int}}
    adjoint_transferlist::Vector{Vector{Int}}
    transferplan::Vector{Vector{PlannedTransfer{Complex{R}}}}
    nodeisfresh::Vector{Bool}
    leafnodeindices::Vector{Int}
    minreceivelevel::Int
    minsourcetranslationlevel::Int
    receivetranslationnodes::Vector{Int}
    firetranslationnodes::Vector{Int}
    aggregationlist::Vector{Vector{Int}}
    disaggregationlist::Vector{Vector{Int}}
    bvector::Vector{Complex{R}}
    verbose::Bool
    tmpmatrix::Matrix{Complex{R}}
end