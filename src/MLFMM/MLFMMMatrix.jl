include("beastglue.jl")
include("MLFMMTree.jl")
# TODD: docstrings
struct MLFMMSource{M<:ResampleMap,Y<:SphereSamplingStrategy,C,X<:MLFMMTree} <:
       AntennaFieldRepresentation{Radiated,C}
    tree::X
    expectedaccuracy::Real
    wavenumber::Real
    nodefarfields::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    basisfunctionfarfields::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    levelcutoffparameters::Vector{Int}
    levelresamplemaps::Vector{M}
    phaseshifttoparent::Array{Matrix{C},2}
    globalphaseshift::Matrix{C}
    nodeisfresh::Vector{Bool}
    nodeisoccupied::Vector{Bool}
    leafnodeindices::Vector{Int}
    aggregationlist::Vector{Vector{Int}}
    rootnode::Int
    buffer::Vector{C}
    outputbuffer::Vector{C}
    verbose::Bool
end
function asvector(p::MLFMMSource)
    return p.buffer
end
Base.size(p::MLFMMSource) = size(asvector(p))
Base.getindex(p::MLFMMSource, i) = Base.getindex(asvector(p), i)
Base.setindex!(p::MLFMMSource, i, v) = Base.setindex!(asvector(p), i, v)
function Base.similar(p::MLFMMSource)
    return deepcopy(p)
end

function MLFMMSource(
    basisfunctions, #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions) 
    wavenumber::T;
    expectedaccuracy=T(1e-3),
    verbose=false,
    minhalfsize=π / (2 * wavenumber),
    orderθ=8,
    orderϕ=8,
    samplingtype::Type{S}=GaussLegendreθRegularϕSampling,
) where {T<:Real,S<:SphereSamplingStrategy}
    verbose &&
        @info "----------------------\n   Assemble MLFMM source \n----------------------------"

    verbose && @info "Initialize source tree"
    points = _getpoints(basisfunctions)
    tree = _initialize_tree(points, minhalfsize)
    levelcutoffparameters =
        _initializelevelcutoffparameters(tree, expectedaccuracy, wavenumber)
    nodeisfresh = [false for _ = 1:length(tree.nodes)]

    verbose && @info "Allocate  node buffers"
    nodeisoccupied = _findoccupiednodes(tree)
    rootnode = _findrootnode(tree, nodeisoccupied)

    min_aggregationlevel = level(tree, rootnode)

    nodefarfields = _allocatenodepattern(
        PlaneWaveExpansion{Radiated,samplingtype,Complex{T}},
        tree,
        levelcutoffparameters,
        nodeisoccupied,
        wavenumber,
    )
    sourceoutputbuffer = nodefarfields[rootnode].buffer

    verbose && @info "Assemble leaf patterns"
    leafnodeindices = leafs(tree)
    basisfunctionfarfields = _initializebasisfunctionfarfields(
        tree,
        basisfunctions,
        levelcutoffparameters[end],
        wavenumber;
        verbose=verbose,
        samplingtype=samplingtype,
    )
    inputbuffer = Vector{Complex{T}}(undef, length(basisfunctions))
    inputbuffer .= basisfunctions

    verbose && @info "Assemble resample maps"
    levelresamplemaps = _initializelevelresamplemaps(
        T,
        tree,
        orderθ,
        orderϕ,
        levelcutoffparameters;
        samplingtype=samplingtype,
    )
    phaseshifttoparent = _initializephaseshifttoparent(
        tree,
        levelcutoffparameters,
        T(wavenumber);
        samplingtype=samplingtype,
    )
    L = levelcutoffparameters[min_aggregationlevel]

    sampling = _standardsampling(samplingtype, L)
    θs, ϕs = samples(sampling)
    nθ, nϕ = length(θs), length(ϕs)
    R = center(tree, rootnode)
    globalphaseshift = Matrix{Complex{T}}(undef, nθ, nϕ)
    _phaseshiftmatrix!(globalphaseshift, -R, wavenumber, sampling)

    numlevels = length(levels(tree))
    # aggregationlist = Vector{Vector{Int}}(undef, numlevels)
    aggregationlist = [Vector{Int}([]) for _ = 1:numlevels]

    verbose && println("------------------------------")

    return MLFMMSource{eltype(levelresamplemaps),samplingtype,Complex{T},typeof(tree)}(
        tree,
        expectedaccuracy,
        wavenumber,
        nodefarfields,
        basisfunctionfarfields,
        levelcutoffparameters,
        levelresamplemaps,
        phaseshifttoparent,
        globalphaseshift,
        nodeisfresh,
        nodeisoccupied,
        leafnodeindices,
        aggregationlist,
        rootnode,
        inputbuffer,
        sourceoutputbuffer,
        verbose,
    )

end

# TODD: docstrings
struct MLFMMReceive{M<:ResampleMap,Y<:SphereSamplingStrategy,C,X<:MLFMMTree}
    tree::X
    expectedaccuracy::Real
    wavenumber::Real
    nodespectra::Vector{PlaneWaveExpansion{Incident,Y,C}}
    testfunctionfarfields::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    levelcutoffparameters::Vector{Int}
    levelresamplemaps::Vector{M}
    phaseshifttoparent::Array{Matrix{C},2}
    nodeisfresh::Vector{Bool}
    nodeisoccupied::Vector{Bool}
    disaggregationlist::Vector{Vector{Int}}
    rootnode::Int
    leafnodeindices::Vector{Int}
    verbose::Bool
    buffer::Vector{C}
end


# TODD: docstrings
struct MLFMMTransmitMap{
    A<:AntennaFieldRepresentation,
    F<:FieldSampling,
    M<:ResampleMap,
    Y<:SphereSamplingStrategy,
    C<:Complex,
    X<:MLFMMTree,
    TP<:AbstractTransfer,
} <: TransmitMap{A,F,C}
    inputbuffer::Vector{C}
    outputbuffer::Vector{C}
    sourcestruct::MLFMMSource{M,Y,C,X}
    receivestruct::MLFMMReceive{M,Y,C,X}
    transferlist::Vector{Vector{Int}}
    adjoint_transferlist::Vector{Vector{Int}}
    transferplan::Vector{Vector{TP}}
    firetranslationnodes::Vector{Int}
    mintranslationlevel::Int
    receivetranslationnodes::Vector{Int}
    verbose::Bool
    # tmpmatrix::Matrix{C}
end

Base.size(p::MLFMMTransmitMap) = (length(p.outputbuffer), length(p.inputbuffer))

function MLFMMTransmitMap(
    basisfunctions, #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions)
    fieldsampling::IrregularFieldSampling,
    wavenumber::T;
    expectedaccuracy=T(1e-3),
    verbose=false,
    minhalfsize=π / (2 * wavenumber),
    orderθ=8,
    orderϕ=8,
    samplingtype::Type{S}=GaussLegendreθRegularϕSampling,
    num_bufferboxes::Integer=1,
    transfertype=PlannedTransfer{samplingtype,Complex{T},T},
) where {T<:Real,S<:SphereSamplingStrategy}
    C = Complex{T}

    verbose && @info "Initialize  MLFMM tree"
    sourcepoints = _getpoints(basisfunctions)
    receivepoints = _getpoints(fieldsampling)
    points = [sourcepoints; receivepoints]
    sourcetree = _initialize_tree(points, minhalfsize)
    leafnodeindices = leafs(sourcetree)
    receivetree = deepcopy(sourcetree)
    for leafnode in leafnodeindices
        # Threads.@threads for leafnode in A.leafnodeindices
        sourceindices =
            findall(x -> x <= length(sourcepoints), sourcetree(leafnode).data.values)
        receiveindices =
            findall(x -> x > length(sourcepoints), receivetree(leafnode).data.values)

        deleteat!(sourcetree(leafnode).data.values, receiveindices)
        deleteat!(receivetree(leafnode).data.values, sourceindices)
        receivetree(leafnode).data.values .-= length(sourcepoints)

    end


    levelcutoffparameters =
        _initializelevelcutoffparameters(sourcetree, expectedaccuracy, wavenumber)

    transmitnodeisoccupied = _findoccupiednodes(sourcetree)
    receivenodeisoccupied = _findoccupiednodes(receivetree)

    sourcerootnode = _findrootnode(sourcetree, transmitnodeisoccupied)
    receiverootnode = _findrootnode(receivetree, receivenodeisoccupied)

    transmitnodeisfresh = [false for k = 1:length(sourcetree.nodes)]
    receivenodeisfresh = [false for k = 1:length(receivetree.nodes)]

    message = string(
        " ",
        length(levels(sourcetree)),
        " levels,  ",
        numberofnodes(sourcetree),
        " boxes",
    )
    verbose && @info message
    message = string(sum(transmitnodeisoccupied), " boxes contain sources")
    verbose && @info message
    message = string(
        "Box ",
        sourcerootnode,
        " at level ",
        level(sourcetree, sourcerootnode),
        " contains all sources",
    )
    verbose && @info message
    message = string(sum(receivenodeisoccupied), " boxes contain probes")
    verbose && @info message
    message = string(
        "Box ",
        receiverootnode,
        " at level ",
        level(receivetree, receiverootnode),
        " contains all probes",
    )
    verbose && @info message


    verbose && println()
    verbose && @info "Allocate node patterns"

    nodefarfields = _allocatenodepattern(
        PlaneWaveExpansion{Radiated,samplingtype,Complex{T}},
        sourcetree,
        levelcutoffparameters,
        transmitnodeisoccupied,
        wavenumber,
    )
    message = string(Base.format_bytes(Base.summarysize(nodefarfields)), " in memory")
    verbose && @info message

    verbose && println()
    verbose && @info "Allocate node  spectra"

    nodespectra = _allocatenodepattern(
        PlaneWaveExpansion{Incident,samplingtype,Complex{T}},
        receivetree,
        levelcutoffparameters,
        receivenodeisoccupied,
        wavenumber,
    )
    message = string(Base.format_bytes(Base.summarysize(nodespectra)), " in memory")
    verbose && @info message

    verbose && println()
    verbose && @info "Compute basis patterns"

    L = levelcutoffparameters[end]
    samplingstrategy = _standardsampling(S, L)

    basisfunctionfarfields = _initializebasisfunctionfarfields(
        sourcetree,
        basisfunctions,
        samplingstrategy,
        wavenumber;
        verbose=verbose,
    )
    inputbuffer = Vector{Complex{T}}(undef, length(basisfunctions))
    inputbuffer .= basisfunctions
    message = string(length(inputbuffer), " basis functions")
    verbose && @info message
    message =
        string(Base.format_bytes(Base.summarysize(basisfunctionfarfields)), " in memory")
    verbose && @info message

    outputbuffer = deepcopy(asvector(fieldsampling))
    verbose && println()
    verbose && @info "Compute probe patterns"
    message = string(length(fieldsampling), " sampling positions")
    verbose && @info message

    testfunctionfarfields = _initializebasisfunctionweightingpatterns(
        receivetree,
        fieldsampling,
        samplingstrategy,
        verbose,
    )
    message =
        string(Base.format_bytes(Base.summarysize(testfunctionfarfields)), " in memory")
    verbose && @info message

    verbose && println()
    verbose && @info "Assemble resample maps"
    levelresamplemaps = _initializelevelresamplemaps(
        T,
        sourcetree,
        orderθ,
        orderϕ,
        levelcutoffparameters;
        samplingtype=samplingtype,
    )
    phaseshifttoparent = _initializephaseshifttoparent(
        sourcetree,
        levelcutoffparameters,
        T(wavenumber);
        samplingtype=samplingtype,
    )
    sumsize = Base.summarysize(levelresamplemaps) + Base.summarysize(phaseshifttoparent)
    message = string(Base.format_bytes(sumsize), " in memory")
    verbose && @info message

    verbose && println()
    verbose && @info "Initialize   transfers"

    transferlist,
    adjoint_transferlist,
    mintranslationlevel,
    receivetranslationnodes,
    firetranslationnodes,
    transferplan,
    uniquetransfers = _initialize_transfers(
        sourcetree,
        receivetree,
        num_bufferboxes,
        transmitnodeisoccupied,
        receivenodeisoccupied,
        nodefarfields,
        nodespectra,
        levelcutoffparameters;
        transfertype=transfertype,
    )

    numtransfers = 0
    for transfers in transferlist
        numtransfers += length(transfers)
    end

    message = string(" ", numtransfers, " planned transfers")
    verbose && @info message

    for k in levels(receivetree)
        message = string(length(uniquetransfers[k]), " unique transfers on level ", k)
        verbose && @info message
    end
    message = string(Base.format_bytes(Base.summarysize(uniquetransfers)), " in memory")
    verbose && @info message

    aggregationlist =
        _initialize_aggregationlist(sourcetree, transferlist, transmitnodeisoccupied)
    disaggregationlist =
        _initialize_disaggregationlist(receivetree, transferlist, receivenodeisoccupied)

    sampling = nodefarfields[sourcerootnode].samplingstrategy
    θs, ϕs = samples(sampling)
    nθ, nϕ = length(θs), length(ϕs)
    R = center(sourcetree, sourcerootnode)
    globalphaseshift = Matrix{Complex{T}}(undef, nθ, nϕ)
    _phaseshiftmatrix!(globalphaseshift, -R, wavenumber, sampling)

    verbose && println("------------------------------")

    sourcestruct =
        MLFMMSource{typeof(levelresamplemaps[end]),samplingtype,C,typeof(sourcetree)}(
            sourcetree,
            expectedaccuracy,
            wavenumber,
            nodefarfields,
            basisfunctionfarfields,
            levelcutoffparameters,
            levelresamplemaps,
            phaseshifttoparent,
            globalphaseshift,
            transmitnodeisfresh,
            transmitnodeisoccupied,
            leafs(sourcetree),
            aggregationlist,
            sourcerootnode,
            inputbuffer,
            nodefarfields[sourcerootnode].buffer,
            verbose,
        )

    receivestruct =
        MLFMMReceive{typeof(levelresamplemaps[end]),samplingtype,C,typeof(sourcetree)}(
            receivetree,
            expectedaccuracy,
            wavenumber,
            nodespectra,
            vec(testfunctionfarfields),
            levelcutoffparameters,
            levelresamplemaps,
            phaseshifttoparent,
            receivenodeisfresh,
            receivenodeisoccupied,
            disaggregationlist,
            receiverootnode,
            leafnodeindices,
            verbose,
            outputbuffer
        )

    transmap = MLFMMTransmitMap{
        typeof(basisfunctions),
        typeof(fieldsampling),
        typeof(levelresamplemaps[end]),
        samplingtype,
        C,
        typeof(sourcetree),
        transfertype,
    }(
        inputbuffer,
        outputbuffer,
        sourcestruct,
        receivestruct,
        transferlist,
        adjoint_transferlist,
        transferplan,
        firetranslationnodes,
        mintranslationlevel,
        receivetranslationnodes,
        verbose,
        # tmpmatrix::Matrix{C}
    )

    return transmap

end

"""
    _findrootnode(tree, nodeisoccupied)

    Returns node on highest level (i.e. smallest box) which contains all leaf elements
"""
function _findrootnode(tree, nodeisoccupied)
    occupiednodesatlevel = [Int[] for _ in levels(tree)]
    rootnode = 0
    for level in levels(tree)
        for node in nodesatlevel(tree, level)
            nodeisoccupied[node] && append!(occupiednodesatlevel[level], node)
        end
        if length(occupiednodesatlevel[level]) == 1
            rootnode = maximum([rootnode, Int(level)])
        end
    end

    return rootnode
end

"""
    _initialize_aggregationlist(sourcetree:MLFMMTree, transmitnodeisoccupied)

Return a list of nodes per level which shall be aggregated to be able to perform all transfers provided in `transferlist`.
"""
function _initialize_aggregationlist(sourcetree, transferlist, transmitnodeisoccupied)
    # sourcetree = MLFMMTrees.tree(TXtree)
    numlevels = length(levels(sourcetree))
    # aggregationlist = Vector{Vector{Int}}(undef, numlevels)
    aggregationlist = [Vector{Int}([]) for _ = 1:numlevels]
    for transmitlist in transferlist
        for translatenode in transmitlist
            for childnode in DepthFirstIterator(sourcetree, translatenode)
                isleaf(sourcetree, childnode) && continue
                !(transmitnodeisoccupied[childnode]) && continue
                lvl = level(sourcetree, childnode)
                push!(aggregationlist[lvl], childnode)
            end
        end
    end

    for level in eachindex(aggregationlist)
        aggregationlist[level] = sort!(unique(aggregationlist[level]))
    end
    return aggregationlist
end

"""
    _initialize_disaggregationlist(RXtree::MLFMMTrees.AbstractMLFMMTree, transferlist)

Return a list of nodes per level which shall be disaggregated to be able to perform all  disaggregations from nodes which receive transfers provided in `transferlist`.
"""
function _initialize_disaggregationlist(receivetree, transferlist, receivenodeisoccupied)
    # receivetree = MLFMMTrees.tree(RXtree)
    numlevels = length(levels(receivetree))
    # disaggregationlist = Vector{Vector{Int}}(undef, numlevels)
    disaggregationlist = [Vector{Int}([]) for k = 1:numlevels]
    for receivenode in eachindex(transferlist)
        if transferlist[receivenode] != []
            for childnode in DepthFirstIterator(receivetree, receivenode)
                isleaf(receivetree, childnode) && continue
                !(receivenodeisoccupied[childnode]) && continue
                lvl = level(receivetree, childnode)
                push!(disaggregationlist[lvl], childnode)
            end
        end
    end
    for level in eachindex(disaggregationlist)
        disaggregationlist[level] = sort!(unique(disaggregationlist[level]))
    end
    return disaggregationlist

end


function _initialize_transfers(
    sourcetree,
    receivetree,
    numbufferboxes,
    transmitnodeisoccupied,
    receivenodeisoccupied,
    nodefarfields,
    nodespectra,
    levelcutoffparameters;
    transfertype::P=PlannedTransfer{
        typeof(nodefarfields[1].sampling),
        eltype(nodefarfields[1].EθEϕ),
        typeof(nodefarfields[1].k0),
    },
) where {P<:Type{<:AbstractTransfer}}

    transferlist = [Int[] for _ = 1:numberofnodes(receivetree)]
    adjoint_transferlist = [Int[] for _ = 1:numberofnodes(sourcetree)]
    mintranslationlevel = typemax(Int)
    uniquetransfers = [transfertype[] for _ in levels(receivetree)]

    transferplan = [
        Vector{transfertype}(undef, length(sourcetree.nodes)) for
        _ in eachindex(transferlist)
    ]

    Pℓstorage =
        Vector{typeof(nodefarfields[1].wavenumber)}(undef, levelcutoffparameters[1] + 1)


    for receivenode in DepthFirstIterator(receivetree, root(receivetree))
        !(receivenodeisoccupied[receivenode]) && continue

        receivelevel = level(receivetree, receivenode)

        for sourcenode in nodesatlevel(sourcetree, receivelevel)
            !(transmitnodeisoccupied[sourcenode]) && continue

            if _transfercanhappen(
                sourcenode,
                receivenode,
                sourcetree,
                numbufferboxes=numbufferboxes,
            )
                append!(transferlist[receivenode], sourcenode)
                append!(adjoint_transferlist[sourcenode], receivenode)
                mintranslationlevel = minimum([Int(receivelevel), mintranslationlevel])
                boxhalfsize = halfsize(receivetree, receivenode)
                transvector =
                    center(sourcetree, sourcenode) - center(receivetree, receivenode)
                isnew = true
                for uniquetransfer in uniquetransfers[receivelevel]
                    uniquetransfervector = gettransfervector(uniquetransfer)
                    if norm(transvector - uniquetransfervector) < 1e-3 * boxhalfsize
                        isnew = false

                        transferplan[receivenode][sourcenode] = uniquetransfer
                        break
                    end
                end

                if isnew

                    receivesampling = nodespectra[receivenode].samplingstrategy
                    sourcesampling = nodefarfields[sourcenode].samplingstrategy
                    θs, ϕs = samples(receivesampling) .- samples(sourcesampling)

                    maximum(abs.(θs)) > 1e-13 && DimensionMismatch(
                        "Samplings before and after Transfer don't match.",
                    )
                    maximum(abs.(ϕs)) > 1e-13 && DimensionMismatch(
                        "Samplings before and after Transfer don't match.",
                    )


                    # newtransfer = 
                    # _initialize_plannedtransfer!(
                    #     Pℓstorage,
                    #     transvector, 
                    #     getwavenumber(nodefarfields[sourcenode]), 
                    #     nodefarfields[sourcenode].samplingstrategy,
                    #     levelcutoffparameters[receivelevel])

                    # transferplan[receivenode][sourcenode] = newtransfer 
                    # append!(uniquetransfers[receivelevel], newtransfer)

                    transferplan[receivenode][sourcenode] = _initialize_plannedtransfer!(
                        Pℓstorage,
                        transvector,
                        getwavenumber(nodefarfields[sourcenode]),
                        nodefarfields[sourcenode].samplingstrategy,
                        levelcutoffparameters[receivelevel],
                        multiplyweights=true,
                    )

                    append!(
                        uniquetransfers[receivelevel],
                        [transferplan[receivenode][sourcenode]],
                    )
                end
            end

        end

    end
    receivetranslationnodes = [k for k in eachindex(transferlist) if transferlist[k] != []]
    firetranslationnodes =
        [k for k in eachindex(adjoint_transferlist) if adjoint_transferlist[k] != []]
    return transferlist,
    adjoint_transferlist,
    mintranslationlevel,
    receivetranslationnodes,
    firetranslationnodes,
    transferplan,
    uniquetransfers
end

"""
    _transfercanhappen(sourcenode, receivenode, tree, numbufferboxes=1)
Returns `true` if a transfer can happen between the `sourcenode` and `receivenode` of the `receivetree` and `false` otherwise.

A transfer can happen if `sourcenode` and `receivenode` are far but their parents are near.

Inputs: 
- `sourcenode` : Index of the source node
- `receivenode` : Index of the receive node
- `tree` : Octree structure containing `sourcenode`
- `numbufferboxes` : minimum number of empty boxes between `sourcenode` and `receivenode` to count as far from each other.  

"""
function _transfercanhappen(sourcenode, receivenode, tree; numbufferboxes=1)
    numbufferboxes = maximum([one(typeof(numbufferboxes)), numbufferboxes])
    sourcelevel = level(tree, sourcenode)
    receivelevel = level(tree, receivenode)
    sourcelevel < 3 && return false
    sourcelevel != receivelevel && return false

    isnearmlfmmbox(
        center(tree, sourcenode),
        center(tree, receivenode),
        halfsize(tree, sourcenode),
        numbufferboxes,
    ) && return false

    sourceparent = parent(tree, sourcenode)
    receiveparent = parent(tree, receivenode)

    isfarmlfmmbox(
        center(tree, sourceparent),
        center(tree, receiveparent),
        halfsize(tree, sourceparent),
        numbufferboxes,
    ) && return false

    return true

end

"""
    isnearmlfmmbox(
    center_a,
    center_b,
    halfsize,
    numbufferboxes,
)

Returns 'true' if each any coordinate of `Rvec = (center_a - center_b)` is smaller than `2 * (bufferboxes + 1 ) * halfsize` and 'false' otherwise.

# Inputs:
- `center_a`: Coordinate vector of the center of the first box
- `center_b`: Coordinate vector of the center of the second box
- `halfsize`: Half of the sidelength of the boxes (sidelength for first and second box are assumed to be equal)
- 'numbufferboxes': minimum number of boxes which must lie between first and second box such that the boxes count as "far" 
"""
function isnearmlfmmbox(center_a, center_b, halfsize, numbufferboxes)
    Rvec = (center_a - center_b) ./ (2 * halfsize)

    threshold = (numbufferboxes + 0.5)^2

    for coordinate in Rvec .^ 2
        if coordinate > threshold
            return false
        end
    end
    return true
end


"""
    isfarmlfmmbox(
    center_a,
    center_b,
    halfsize,
    numbufferboxes,
)

Returns 'false' if each any coordinate of `Rvec = (center_a - center_b)` is smaller than `2 * (bufferboxes + 1 ) * halfsize` and 'true' otherwise.

# Inputs:
- `center_a`: Coordinate vector of the center of the first box
- `center_b`: Coordinate vector of the center of the second box
- `halfsize`: Half of the sidelength of the boxes (sidelength for first and second box are assumed to be equal)
- 'numbufferboxes': minimum number of boxes which must lie between first and second box such that the boxes count as "far" 
"""
function isfarmlfmmbox(center_a, center_b, halfsize, numbufferboxes)
    return !(isnearmlfmmbox(center_a, center_b, halfsize, numbufferboxes))
end

"""
    _findoccupiednodes(tree; minlevel =1)

Return a voctor of Boolean values to indicate which nodes of the tree contain values.
"""
function _findoccupiednodes(tree; minlevel=1)
    nodeisoccupied = [false for k = 1:length(tree.nodes)]

    levels = AntennaFieldRepresentations.levels(tree)
    for level in reverse(maximum([minlevel, 1]):length(levels))
        for parentnode::Int in nodesatlevel(tree, level)
            if isleaf(tree, parentnode)
                if !(isempty(tree(parentnode).data.values))
                    nodeisoccupied[parentnode] = true
                end
            end
            for child in children(tree, parentnode)
                if nodeisoccupied[child]
                    nodeisoccupied[parentnode] = true
                end
            end

        end
    end
    return nodeisoccupied
end

"""
    _getpoints(basisfunctions)

Return list of `SVector{3}` representing the reference (i.e., center) locations of the `basisfunctions`
"""
function _getpoints(basisfunctions::SurfaceCurrentDensity)
    return BEAST.positions(basisfunctions.functionspace)
end
function _getpoints(functionspace)
    return SVector{3}.(functionspace.points)
end
function _getpoints(functionspace::BEAST.Space{R}) where {R<:Real}
    return BEAST.positions(functionspace)
end
function _getpoints(functionspace::DipoleArray)
    return SVector{3}.(functionspace.positions)
end
function _getpoints(sampling::IrregularFieldSampling)
    return SVector{3}.(vec(sampling.positions))
end

"""
    _initialize_tree(points::AbstractArray{SVector{3,R}}, minhalfsize::R) where{R<:Real}

Assemble 'MLFMMTree' based on point cloud and half length of leave boxes
"""
function _initialize_tree(
    points::AbstractArray{SVector{3,R}},
    minhalfsize::R,
) where {R<:Real}
    rootcenter, rootsize = getboundingbox(points)

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2 = log2(rootsize / minhalfsize)
    halfsize_powerof2 = maximum([0, Integer(ceil(halfsize_powerof2))])
    rootsize = 2^halfsize_powerof2 * minhalfsize

    root_center = SVector{3,Float64}(rootcenter)
    return AntennaFieldRepresentations.MLFMMTree(root_center, points, rootsize, minhalfsize)
end

"""
    _cutoffparameter(diameter, expectedaccuracy, k0)

Return cutoff parameter (i.e., maximum spherical mode index L) to represent sources inside sphere of given `diameter` with expected accuracy.
"""
function _cutoffparameter(diameter, expectedaccuracy, k0)
    kd = k0 * diameter
    digitsofaccuracy = -log10(expectedaccuracy)
    L = ceil(Int, kd + 1.8 * digitsofaccuracy^(2 / 3) * (kd)^(1 / 3)) + 2
    return L
end

"""
    _initializelevelcutoffparameters(tree::MLFMMTrees.AbstractMLFMMTree, expectedaccuracy, wavenumber)

Return list of cutoff parameters for each level in `tree` to represent the node farfields with expected accuracy.
"""
function _initializelevelcutoffparameters(tree::MLFMMTree, expectedaccuracy, wavenumber)
    levels = AntennaFieldRepresentations.levels(tree)
    cutoffparameters = Array{Int}(undef, length(levels))

    for level in levels
        node = tree.nodes[nodesatlevel(tree, level)][1]
        halfsize = node.data.halfsize
        diameter = sqrt(3 * halfsize^2) * 2
        cutoffparameters[level] = _cutoffparameter(diameter, expectedaccuracy, wavenumber)
    end

    return cutoffparameters
end

"""
    _allocatenodepattern(W::Type{PWS}, tree::MLFMMTree, cutoffparameters::Vector{<:Integer}, [minlevelel::Int]) where {PWS<:PlaneWaveSpectrum}

Allocate memory for intermediate representation of far field patterns for each node in tree. 
"""
# function _allocatenodepattern(
#     W::Type{PlaneWaveExpansion{P,S,C}},
#     tree::MLFMMTree,
#     cutoffparameters::Vector{<:Integer},
#     wavenumber;
#     minlevel::Integer = 1,
# ) where {P<:PropagationType,S<:SphereSamplingStrategy,C}
#     nodepattern = Vector{W}(undef, length(tree.nodes))

#     levels = AntennaFieldRepresentations.levels(tree)
#     for level = minlevel:length(levels)
#         L = cutoffparameters[level]
#         samplingstrategy = _standardsampling(S, L)
#         θs, ϕs = samples(samplingstrategy)
#         for node in AntennaFieldRepresentations.nodesatlevel(tree, level)
#             Eθϕ = zeros(C, length(θs), length(ϕs), 2)
#             buffer = vec(Eθϕ)
#             nodepattern[node] = W(samplingstrategy, Eθϕ, wavenumber, buffer)
#         end
#     end

#     return nodepattern
# end
function _allocatenodepattern(
    W::Type{PlaneWaveExpansion{P,S,C}},
    tree::MLFMMTree,
    cutoffparameters::Vector{<:Integer},
    nodeisoccupied::Vector{Bool},
    wavenumber;
    minlevel::Integer=1,
) where {P<:PropagationType,S<:SphereSamplingStrategy,C}
    nodepattern = Vector{W}(undef, length(tree.nodes))


    levels = AntennaFieldRepresentations.levels(tree)

    for level = minlevel:length(levels)
        L = cutoffparameters[level]
        samplingstrategy = _standardsampling(S, L)
        θs, ϕs = samples(samplingstrategy)
        for node in AntennaFieldRepresentations.nodesatlevel(tree, level)
            !(nodeisoccupied[node]) && continue
            Eθϕ = zeros(C, length(θs), length(ϕs), 2)
            buffer = vec(Eθϕ)
            nodepattern[node] = W(samplingstrategy, Eθϕ, wavenumber, buffer)
        end
    end

    return nodepattern
end

"""
    _initializebasisfunctionfarfields(tree::MLFMMTrees.AbstractMLFMMTree, basisfunctions::S, cutoffparameters::Vector{I}, k0::R; verbose::Bool=false) where {R<:Real, S<:AbstractSurfaceCurrentDensity, I<:Int}

Return list of far fields for every basis function on leaf level.
"""
function _initializebasisfunctionfarfields(
    tree::MLFMMTree,
    basisfunctions::S,
    cutoffparameter::Integer,
    k0::T;
    samplingtype::Type{Y}=GaussLegendreθRegularϕSampling,
    verbose::Bool=false,
) where {T<:Real,S,Y<:SphereSamplingStrategy}
    L = cutoffparameter
    samplingstrategy = _standardsampling(samplingtype, L)
    return _initializebasisfunctionfarfields(
        tree,
        basisfunctions,
        samplingstrategy,
        k0;
        verbose=verbose,
    )

end
function _initializebasisfunctionfarfields(
    tree::MLFMMTree,
    basisfunctions::S,
    samplingstrategy::Y,
    k0::T;
    verbose::Bool=false,
) where {T<:Real,S,Y<:SphereSamplingStrategy}

    θvec, ϕvec = samples(samplingstrategy)
    nθ, nϕ = length(θvec), length(ϕvec)
    basisfunctionfarfields = individualfarfields(basisfunctions, samplingstrategy)
    phaseshiftmatrix = Array{Complex{T}}(undef, nθ, nϕ)
    for leafnode::Int in leafs(tree)
        _phaseshiftmatrix!(phaseshiftmatrix, center(tree, leafnode), k0, samplingstrategy)
        for functionindex::Int in tree(leafnode).data.values
            _eθ(basisfunctionfarfields[functionindex]) .*= (phaseshiftmatrix)
            _eϕ(basisfunctionfarfields[functionindex]) .*= (phaseshiftmatrix)
        end
    end

    return basisfunctionfarfields
end

"""
    individualfarfields(basisfunctions, sampling::Y) where {Y<:SphereSamplingStrategy}

Returns a vector of farfields, one for each coefficient in `basisfunctions`, sampled according to `sampling`.
"""
function individualfarfields(
    basisfunctions::SurfaceCurrentDensity{P,E,B,C},
    sampling::Y,
) where {Y<:SphereSamplingStrategy,P,E,B,C}
    T = real(C)
    θvec, ϕvec = samples(sampling)

    nθ, nϕ = length(θvec), length(ϕvec)
    cost = cos.(θvec)
    sint = sin.(θvec)

    cosp = cos.(ϕvec)
    sinp = sin.(ϕvec)

    eθ =
        SVector{
            3,
        }.([
            [cost[kθ] * cosp[kϕ], cost[kθ] * sinp[kϕ], -sint[kθ]] for kθ in eachindex(θvec),
            kϕ in eachindex(ϕvec)
        ])
    eϕ = SVector{3}.([[-sinp[kϕ], cosp[kϕ], zero(T)] for kϕ in eachindex(ϕvec)])

    pts = [
        point(cosp[kk] * sint[k], sinp[kk] * sint[k], cost[k]) for k in eachindex(θvec),
        kk in eachindex(ϕvec)
    ]

    farfieldmatrices =
        individualcartesianfarfields(basisfunctions, pts)::Matrix{SVector{3,C}}

    basisfunctionfarfields =
        Vector{PlaneWaveExpansion{Radiated,Y,C}}(undef, numfunctions(basisfunctions))
    ff = Matrix{SVector{3,C}}(undef, nθ, nϕ)
    for functionindex in eachindex(basisfunctionfarfields)
        ff .= reshape(view(farfieldmatrices, :, functionindex), nθ, nϕ)

        basisfunctionfarfields[functionindex] = PlaneWaveExpansion(
            Radiated(),
            sampling,
            Matrix{C}(undef, nθ, nϕ),
            Matrix{C}(undef, nθ, nϕ),
            getwavenumber(basisfunctions),
        )
        # Eθ = _eθ(basisfunctionfarfields[functionindex])
        # Eϕ = _eϕ(basisfunctionfarfields[functionindex])
        for kθ in eachindex(θvec), kϕ in eachindex(ϕvec)


            # Eθ[kθ, kϕ] = C(udot(eθ[kθ, kϕ], ff[kθ, kϕ]))
            # Eϕ[kθ, kϕ] = C(udot(eϕ[kϕ], ff[kθ, kϕ]))

            basisfunctionfarfields[functionindex].EθEϕ[kθ, kϕ, 1] =
                C(udot(eθ[kθ, kϕ], ff[kθ, kϕ]))
            basisfunctionfarfields[functionindex].EθEϕ[kθ, kϕ, 2] =
                C(udot(eϕ[kϕ], ff[kθ, kϕ]))
        end

    end
    return basisfunctionfarfields
end
function individualfarfields(
    dipoles::DipoleArray{Radiated,E,T,C},
    sampling::Y,
) where {E,T,C,Y<:SphereSamplingStrategy}
    θvec, ϕvec = samples(sampling)

    nθ, nϕ = length(θvec), length(ϕvec)
    cost = cos.(θvec)
    sint = sin.(θvec)

    cosp = cos.(ϕvec)
    sinp = sin.(ϕvec)


    eθ =
        SVector{
            3,
        }.([
            [cost[kθ] * cosp[kϕ], cost[kθ] * sinp[kϕ], -sint[kθ]] for kθ in eachindex(θvec),
            kϕ in eachindex(ϕvec)
        ])
    eϕ = SVector{3}.([[-sinp[kϕ], cosp[kϕ], zero(T)] for kϕ in eachindex(ϕvec)])
    eᵣ =
        SVector{
            3,
        }.([
            [sint[kθ] * cosp[kϕ], sint[kθ] * sinp[kϕ], cost[kθ]] for kθ in eachindex(θvec),
            kϕ in eachindex(ϕvec)
        ])

    k₀ = dipoles.wavenumber


    basisfunctionfarfields =
        Vector{PlaneWaveExpansion{Radiated,Y,Complex{T}}}(undef, length(dipoles))
    for functionindex in eachindex(dipoles)
        basisfunctionfarfields[functionindex] = PlaneWaveExpansion(
            Radiated(),
            sampling,
            Matrix{Complex{T}}(undef, nθ, nϕ),
            Matrix{Complex{T}}(undef, nθ, nϕ),
            getwavenumber(dipoles),
        )
    end


    for kθ in eachindex(θvec), kϕ in eachindex(ϕvec)
        E_FF =
            C(0.0, -k₀) * _dipolefarfieldscalingfactor(E()) / (4π) .*
            cis.(k₀ * udot.(Ref(eᵣ[kθ, kϕ]), dipoles.positions))
        for (i, dir) in enumerate(dipoles.orientations)
            Epolθ, Epolϕ = _dipoledarfieldpolarization(eθ[kθ, kϕ], eϕ[kϕ], dir, E())
            _eθ(basisfunctionfarfields[i])[kθ, kϕ] = E_FF[i] * Epolθ
            _eϕ(basisfunctionfarfields[i])[kθ, kϕ] = E_FF[i] * Epolϕ
        end
    end

    # phaseshiftmatrix = Array{Complex{T}}(undef, nθ, nϕ)
    # for i in eachindex(dipoles)
    #     _phaseshiftmatrix!(phaseshiftmatrix, dipoles.positions[i], getwavenumber(dipoles), sampling)
    #     _eθ(basisfunctionfarfields[i]) .*= phaseshiftmatrix
    #     _eϕ(basisfunctionfarfields[i]) .*= phaseshiftmatrix
    # end

    return basisfunctionfarfields
end

"""
    _initializebasisfunctionweightingpatterns(tree::MLFMMTrees.AbstractMLFMMTree, functionspace,
    cutoffparameters::Vector{<:Integer}, k0::R, verbose::Bool) where{R<:Real}

Return list of far fields for every test function on leaf level stored in reverted direction feasible for testing incident plane wave spectra. 
"""
function _initializebasisfunctionweightingpatterns(
    tree::MLFMMTree,
    fieldsampling::IFS,
    samplingstrategy::Y,
    verbose::Bool,
) where {IFS<:IrregularFieldSampling,Y<:SphereSamplingStrategy}
    θvec, ϕvec = samples(samplingstrategy)
    θvec .= pi .- θvec
    ϕvec .= mod2pi.(ϕvec .+ pi)

    nθ, nϕ = length(θvec), length(ϕvec)

    C = eltype(fieldsampling)

    probepatterns =
        Vector{PlaneWaveExpansion{Radiated,Y,C}}(undef, length(fieldsampling.probes))



    for (probeID, probe) in enumerate(fieldsampling.probes)
        probepatterns[probeID] = PlaneWaveExpansion(
            Radiated(),
            samplingstrategy,
            zeros(C, nθ, nϕ),
            zeros(C, nθ, nϕ),
            getwavenumber(probe),
        )

        for (θind, θ) in enumerate(θvec)
            for (ϕind, ϕ) in enumerate(ϕvec)
                _eθ(probepatterns[probeID])[θind, ϕind],
                _eθ(probepatterns[probeID])[θind, ϕind] = farfield(probe.aut_field, (θ, ϕ))
            end
        end
    end

    basisfunctionweightingpatterns =
        Array{PlaneWaveExpansion{Radiated,Y,C}}(undef, size(fieldsampling.probeIDs))

    for k in eachindex(basisfunctionweightingpatterns)
        probepattern = probepatterns[fieldsampling.probeIDs[k]]
        χ, θ, ϕ = fieldsampling.eulerangles[k]
        basisfunctionweightingpatterns[k] =
            rotate(probepattern, χ, θ, ϕ; orderθ=6, orderϕ=6)
    end

    k0 = getwavenumber(probepatterns[1])
    phaseshiftmatrix = Array{C}(undef, nθ, nϕ)
    for leafnode::Int in leafs(tree)
        _phaseshiftmatrix!(phaseshiftmatrix, center(tree, leafnode), -k0, samplingstrategy)
        for functionindex::Int in tree(leafnode).data.values
            _eθ(basisfunctionweightingpatterns[functionindex]) .*= (phaseshiftmatrix)
            _eϕ(basisfunctionweightingpatterns[functionindex]) .*= (phaseshiftmatrix)
        end
    end

    return basisfunctionweightingpatterns
end

"""
    _initializelevelresamplemaps(R::Type{<:Real}, tree::MLFMMTree, orderθ::Integer, orderϕ::Integer,cutoffparameters::Vector{<:Integer}; minlevel::Int=0)

Return list of interpolators to aggregate from each level to the respective parent level.
"""
function _initializelevelresamplemaps(
    T::Type{<:Real},
    tree::MLFMMTree,
    orderθ::Integer,
    orderϕ::Integer,
    cutoffparameters::Vector{<:Integer};
    minlevel::Int=0,
    samplingtype::Type{Y}=GaussLegendreθRegularϕSampling,
) where {Y}
    levels = AntennaFieldRepresentations.levels(tree)

    θmaptype = LocalθResampleMap{samplingtype,samplingtype,orderθ,T}
    ϕmaptype = LocalϕResampleMap{samplingtype,samplingtype,orderϕ,T}

    θϕmaptype = θϕResampleMap{θmaptype,ϕmaptype,samplingtype,samplingtype,orderθ,orderϕ,T}

    levelresamplemaps = Vector{θϕmaptype}(undef, length(levels) - 1)
    for level in levels[1:end-1]
        level < minlevel && continue
        Lold = cutoffparameters[level+1]
        Lnew = cutoffparameters[level]

        oldsampling = _standardsampling(samplingtype, Lold)
        newsampling = _standardsampling(samplingtype, Lnew)


        levelresamplemaps[level] = LocalθLocalϕResampleMap(
            newsampling,
            oldsampling;
            orderθ=orderθ,
            orderϕ=orderϕ,
        )
    end
    return levelresamplemaps
end
#TODO: only initialize required phaseshifts and share between source and receive
function _initializephaseshifttoparent(
    tree::MLFMMTree,
    cutoffparameters::Vector{I},
    k0::T;
    minlevel::Int=0,
    samplingtype::Type{Y}=GaussLegendreθRegularϕSampling,
) where {Y,T<:Real,I<:Integer}
    levels = AntennaFieldRepresentations.levels(tree)
    phaseshifttoparent = Array{Matrix{Complex{T}}}(undef, 8, length(levels))
    for level in levels[2:end]
        level < minlevel && continue

        L = cutoffparameters[level-1]

        sampling = _standardsampling(samplingtype, L)
        θs, ϕs = samples(sampling)
        nθ, nϕ = length(θs), length(ϕs)
        # phaseshiftmatrix=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        node = tree.nodes[nodesatlevel(tree, level)[1]]
        halfsize = node.data.halfsize

        R = [halfsize, halfsize, halfsize]
        phaseshifttoparent[1, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[1, level], R, k0, sampling)

        R = [-halfsize, halfsize, halfsize]
        phaseshifttoparent[2, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[2, level], R, k0, sampling)

        R = [halfsize, -halfsize, halfsize]
        phaseshifttoparent[3, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[3, level], R, k0, sampling)

        R = [-halfsize, -halfsize, halfsize]
        phaseshifttoparent[4, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[4, level], R, k0, sampling)

        R = [halfsize, halfsize, -halfsize]
        phaseshifttoparent[5, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[5, level], R, k0, sampling)

        R = [-halfsize, halfsize, -halfsize]
        phaseshifttoparent[6, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[6, level], R, k0, sampling)

        R = [halfsize, -halfsize, -halfsize]
        phaseshifttoparent[7, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[7, level], R, k0, sampling)

        R = [-halfsize, -halfsize, -halfsize]
        phaseshifttoparent[8, level] = Matrix{Complex{T}}(undef, nθ, nϕ)
        _phaseshiftmatrix!(phaseshifttoparent[8, level], R, k0, sampling)
    end
    return phaseshifttoparent
end


"""
    _aggregate_leafnodes!(A::MLFMMSource)

Fill storage for leaf node patterns due to excitation vector `A.xvector`.
"""
function _aggregate_leafnodes!(A::MLFMMSource)
    tree = A.tree

    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        # Threads.@threads for leafnode in A.leafnodeindices
        reset = true
        for functionindex::Int in tree(leafnode).data.values::Vector{Int}
            A.nodefarfields[leafnode].buffer .= _muladd_or_mulreset!(
                A.nodefarfields[leafnode],
                A.basisfunctionfarfields[functionindex],
                A.buffer[functionindex],
                reset=reset,
            )
            reset = false
        end
    end
end
"""
    _adjoint_aggregate_leafnodes!(A::MLFMMSource)

Perform the adjoint operation (i.e., conjugate of transposed operation) to `_aggregate_leafnodes!`
"""
function _adjoint_aggregate_leafnodes!(A::MLFMMSource)
    tree = A.tree
    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        pws = A.nodefarfields[leafnode]
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff = A.basisfunctionfarfields[basisfunctionindex]
            A.buffer[basisfunctionindex] = dot(ff, pws)
        end
    end
end
"""
    _transpose_aggregate_leafnodes!(A::MLFMMSource)

Perform the transpose operation to `_aggregate_leafnodes!`
"""
function _transpose_aggregate_leafnodes!(A::MLFMMSource)
    tree = A.tree
    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        pws = A.nodefarfields[leafnode]
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff = A.basisfunctionfarfields[basisfunctionindex]
            A.buffer[basisfunctionindex] = udot(ff, pws)
        end
    end
end

#TODO: store children per node to remove dependency on ClusteTrees.children 
"""
    _aggregate_children!(A::MLFMMSource, parentnode)

Fill pattern storage of parentnode with aggregated pattern from all its children.
"""
function _aggregate_children!(A::MLFMMSource, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)

    resamplemap = A.levelresamplemaps[level]
    reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        resamplemap.outputbuffer .=
            mul!(resamplemap.outputbuffer, resamplemap, A.nodefarfields[child])
        # _eθ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(
        #     _eθ(A.nodefarfields[parentnode]),
        #     view(resamplemap.outputbuffermat, :, :, 1),
        #     A.phaseshifttoparent[sector, level+1],
        #     reset = reset,
        # )
        # _eϕ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(
        #     _eϕ(A.nodefarfields[parentnode]),
        #     view(resamplemap.outputbuffermat, :, :, 2),
        #     A.phaseshifttoparent[sector, level+1],
        #     reset = reset,
        # )
        _muladd_or_mulreset!(
            _eθ(A.nodefarfields[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            A.phaseshifttoparent[sector, level+1],
            reset=reset,
        )
        _muladd_or_mulreset!(
            _eϕ(A.nodefarfields[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            A.phaseshifttoparent[sector, level+1],
            reset=reset,
        )

        reset = false
    end
end
"""
    _transpose_aggregate_children!(A::MLFMMSource, parentnode)

Perform transposed operation to `_aggregate_children!`
"""
function _transpose_aggregate_children!(A::MLFMMSource, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)

    resamplemap = A.levelresamplemaps[level]
    transpose_resamplemat = transpose(resamplemap)
    # reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        # A.nodefarfields[child] .=
        mul!(A.nodefarfields[child], transpose_resamplemat, resamplemap.outputbuffer)


        # resamplemap.outputbuffer .= mul!(resamplemap.outputbuffer, resamplemap, A.nodefarfields[child])
        # _eθ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eθ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 1),  A.phaseshifttoparent[sector, level+1], reset = reset)
        # _eϕ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eϕ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 2),  A.phaseshifttoparent[sector, level+1], reset = reset)

        # reset = false
    end
end
"""
    _adjoint_aggregate_children!(A::MLFMMSource, parentnode)

Perform adjoint operation to `_aggregate_children!`
"""
function _adjoint_aggregate_children!(A::MLFMMSource, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)

    resamplemap = A.levelresamplemaps[level]
    adjoint_resamplemat = adjoint(resamplemap)
    # reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(A.nodefarfields[parentnode]) .* conj.(A.phaseshifttoparent[sector, level+1])
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(A.nodefarfields[parentnode]) .* conj.(A.phaseshifttoparent[sector, level+1])
        A.nodefarfields[child] .= mul!(
            A.nodefarfields[child].buffer,
            adjoint_resamplemat,
            resamplemap.outputbuffer,
        )


        # resamplemap.outputbuffer .= mul!(resamplemap.outputbuffer, resamplemap, A.nodefarfields[child])
        # _eθ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eθ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 1),  A.phaseshifttoparent[sector, level+1], reset = reset)
        # _eϕ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eϕ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 2),  A.phaseshifttoparent[sector, level+1], reset = reset)

        # reset = false
    end
end

"""
    _aggregate_to_aggregationlist!(A::MLFMMSource)

Aggregate `A` according to `A.aggregationlist`.

`A.aggregationlist` specifies for each level which nodes shall be aggregated.
"""
function _aggregate_to_aggregationlist!(A::MLFMMSource)
    A.verbose && @info "Aggregate node far fields"

    _aggregate_leafnodes!(A)

    aggregationlist = A.aggregationlist

    for level in reverse(eachindex(aggregationlist))
        for node in aggregationlist[level]
            # Threads.@threads for node in aggregationlist[level]
            _aggregate_children!(A, node)
        end
    end

end

"""
    _adjoint_aggregate_to_aggregationlist!(A::MLFMMSource)

Perform adjoint operation (i.e., complex conjugate of transposed operation) of ` _transpose_aggregate_to_aggregationlist!`
"""
function _adjoint_aggregate_to_aggregationlist!(A::MLFMMSource)
    A.verbose && @info "Adjoint aggregate node far fields"

    aggregationlist = A.aggregationlist

    for level in eachindex(aggregationlist)
        aggregationlist[level] == [] && continue
        for sector in 1:8
            # Threads.@threads for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end

        for node in aggregationlist[level]
            # Threads.@threads for node in aggregationlist[level]
            _transpose_aggregate_children!(A, node)
        end

        for sector in 1:8
            # Threads.@threads for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end
    end

    _adjoint_aggregate_leafnodes!(A)

end

"""
    _transpose_aggregate_to_aggregationlist!(A::MLFMMSource)

Perform transpose operation of `_transpose_aggregate_to_aggregationlist!`
"""
function _transpose_aggregate_to_aggregationlist!(A::MLFMMSource)
    A.verbose && @info "Adjoint aggregate node far fields"

    aggregationlist = A.aggregationlist

    for level in eachindex(aggregationlist)
        aggregationlist[level] == [] && continue

        for node in aggregationlist[level]
            # Threads.@threads for node in aggregationlist[level]
            _transpose_aggregate_children!(A, node)
        end
    end

    _transpose_aggregate_leafnodes!(A)

end



"""
    _aggregate_to_minlevel!(A::MLFMMSource, [x::AbstractVector]; min_aggregationlevel::Integer=0)

Aggregate `A` up to min_aggregationlevel. 
"""
function _aggregate_to_minlevel!(A::MLFMMSource, x; min_aggregationlevel::Integer=0)
    A.buffer .= x
    _aggregate_to_minlevel!(A, min_aggregationlevel=min_aggregationlevel)
end
function _aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    _aggregate_leafnodes!(A)

    levels = AntennaFieldRepresentations.levels(tree)
    for level in
        reverse(maximum([min_aggregationlevel, level(tree, A.rootnode), 1]):length(levels))
        for parentnode::Int in nodesatlevel(tree, level)
            !(A.nodeisoccupied[parentnode]) && continue
            isleaf(tree, parentnode) && continue
            _aggregate_children!(A, parentnode)
        end
    end

end
"""
    _adjoint_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)


Perform adjoint operation to `aggregate_to_minlevel!`
"""
function _adjoint_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    levels = AntennaFieldRepresentations.levels(tree)
    for level = maximum([min_aggregationlevel, level(tree, A.rootnode), 1]):length(levels)
        # for parentnode::Int in nodesatlevel(tree, level)
        for parentnode::Int in nodesatlevel(tree, level)
            !(A.nodeisoccupied[parentnode]) && continue
            isleaf(tree, parentnode) && continue
            _adjoint_aggregate_children!(A, parentnode)
        end
    end
    _adjoint_aggregate_leafnodes!(A)

end
"""
    _transpose_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)


Perform transpose operation to `aggregate_to_minlevel!`
"""
function _transpose_aggregate_to_minlevel!(
    A::MLFMMSource;
    min_aggregationlevel::Integer=0,
)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    levels = AntennaFieldRepresentations.levels(tree)
    for level = maximum([min_aggregationlevel, level(tree, A.rootnode), 1]):length(levels)
        # for parentnode::Int in nodesatlevel(tree, level)
        for parentnode::Int in nodesatlevel(tree, level)
            !(A.nodeisoccupied[parentnode]) && continue
            isleaf(tree, parentnode) && continue
            _transpose_aggregate_children!(A, parentnode)
        end
    end
    _transpose_aggregate_leafnodes!(A)
end

"""
    _aggregate_to_farfield!(A::MLFMMSource, [x::AbstractVector]; min_aggregationlevel::Integer=0)

Aggregate `A` to farfield. 
"""
function _aggregate_to_farfield!(A::MLFMMSource, x)
    A.buffer .= x
    _aggregate_to_farfield!(A)
end
function _aggregate_to_farfield!(A::MLFMMSource)
    _aggregate_to_minlevel!(A)
    _eθ(A.nodefarfields[A.rootnode]) .=
        _eθ(A.nodefarfields[A.rootnode]) .* A.globalphaseshift
    _eϕ(A.nodefarfields[A.rootnode]) .=
        _eϕ(A.nodefarfields[A.rootnode]) .* A.globalphaseshift
end
"""
    _adjoint_aggregate_to_farfield!(A::MLFMMSource; min_aggregationlevel::Integer=0)

Perform adjoint operation to `_aggregate_to_farfield!`
"""
function _adjoint_aggregate_to_farfield!(A::MLFMMSource)
    _eθ(A.nodefarfields[A.rootnode]) .=
        _eθ(A.nodefarfields[A.rootnode]) .* conj.(A.globalphaseshift)
    _eϕ(A.nodefarfields[A.rootnode]) .=
        _eϕ(A.nodefarfields[A.rootnode]) .* conj.(A.globalphaseshift)
    _adjoint_aggregate_to_minlevel!(A)
end
"""
    _transpose_aggregate_to_farfield!(A::MLFMMSource; min_aggregationlevel::Integer=0)

Perform transpose operation to `_aggregate_to_farfield!`
"""
function _transpose_aggregate_to_farfield!(A::MLFMMSource)
    _eθ(A.nodefarfields[A.rootnode]) .=
        _eθ(A.nodefarfields[A.rootnode]) .* A.globalphaseshift
    _eϕ(A.nodefarfields[A.rootnode]) .=
        _eϕ(A.nodefarfields[A.rootnode]) .* A.globalphaseshift
    _transpose_aggregate_to_minlevel!(A)
end

function equivalentorder(A::MLFMMSource)
    return A.levelcutoffparameters[A.rootnode]
end


"""
    _disaggregate_leafnodes!(receivestruct::MLFMMReceive)

Test leaf node patterns with all test functions of corresponding leafnode and store result in `receivestruct.bvector`
"""
function _disaggregate_leafnodes!(receivestruct::MLFMMReceive)
    receivetree = receivestruct.tree

    # Threads.@threads for leafnode in receivestruct.leafnodeindices
    for leafnode in receivestruct.leafnodeindices
        !(receivestruct.nodeisoccupied[leafnode]) && continue

        pws = receivestruct.nodespectra[leafnode]

        for probeindex::Int in receivetree(leafnode).data.values::Vector{Int}
            ff = receivestruct.testfunctionfarfields[probeindex]

            receivestruct.buffer[probeindex] = udot(pws, ff)
        end
    end

end


"""
    _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)

Perform ajoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate_leafnodes!`
"""
function _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)
    receivetree = A.tree

    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        reset = true
        for functionindex in receivetree(leafnode).data.values
            if !reset
                _eθ(A.nodespectra[leafnode]) .+=
                    A.buffer[functionindex] .*
                    conj.(_eθ(A.testfunctionfarfields[functionindex]))
                _eϕ(A.nodespectra[leafnode]) .+=
                    A.buffer[functionindex] .*
                    conj.(_eϕ(A.testfunctionfarfields[functionindex]))
            else
                _eθ(A.nodespectra[leafnode]) .=
                    A.buffer[functionindex] .*
                    conj.(_eθ(A.testfunctionfarfields[functionindex]))
                _eϕ(A.nodespectra[leafnode]) .=
                    A.buffer[functionindex] .*
                    conj.(_eϕ(A.testfunctionfarfields[functionindex]))
                reset = false
            end
        end
    end

end

"""
    _transpose_disaggregate_leafnodes!(A::MLFMMReceive)

Perform transposed operation to `disaggregate_leafnodes!`
"""
function _transpose_disaggregate_leafnodes!(A::MLFMMReceive)
    receivetree = A.tree

    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        reset = true
        for functionindex::Int in receivetree(leafnode).data.values::Vector{Int}
            A.nodespectra[leafnode].buffer .= _muladd_or_mulreset!(
                A.nodespectra[leafnode],
                A.testfunctionfarfields[functionindex],
                A.buffer[functionindex],
                reset=reset,
            )
            reset = false
        end
    end

end

"""
    _disaggregate_children!(receivestruct, parentnode)

Disaggregate pattern of parentnode to all its children and store the resulting patterns in `receivestruct.nodespectra[child]`
"""
function _disaggregate_children!(receivestruct::MLFMMReceive, parentnode)

    receivetree = receivestruct.tree


    lvl = AntennaFieldRepresentations.level(receivetree, parentnode)

    resamplemap = receivestruct.levelresamplemaps[lvl]
    transpose_resamplemat = transpose(resamplemap)
    # reset = true

    for child in children(receivetree, parentnode)
        !(receivestruct.nodeisoccupied[child]) && continue

        sector = receivetree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(receivestruct.nodespectra[parentnode]) .* receivestruct.phaseshifttoparent[sector, lvl+1]
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(receivestruct.nodespectra[parentnode]) .* receivestruct.phaseshifttoparent[sector, lvl+1]
        # receivestruct.nodespectra[child] .=
        mul!(receivestruct.nodespectra[child], transpose_resamplemat, resamplemap.outputbuffer)
        # receivestruct.nodeisfresh[child] = true
    end
end

"""
    _transpose_disaggregate_children!(A, parentnode)

Perform transposed operation to `disaggrgate_children!`
"""
function _transpose_disaggregate_children!(A, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)

    resamplemap = A.levelresamplemaps[level]
    reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        resamplemap.outputbuffer .=
            mul!(resamplemap.outputbuffer, resamplemap, A.nodespectra[child])

        _muladd_or_mulreset!(
            _eθ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            A.phaseshifttoparent[sector, level+1],
            reset=reset,
        )
        _muladd_or_mulreset!(
            _eϕ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            A.phaseshifttoparent[sector, level+1],
            reset=reset,
        )

        reset = false
    end
end

"""
    _adjoint_disaggregate_children!(A, parentnode)

Perform adjoint operation to `disaggrgate_children!`
"""
function _adjoint_disaggregate_children!(A, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)


    # resamplemap = conj.(A.levelresamplemaps[level]) -> since resamplemap is real, conj.(conj.(A.levelresamplemaps[level]) == A.levelresamplemaps[level])
    resamplemap = A.levelresamplemaps[level]
    reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        resamplemap.outputbuffer .=
            mul!(resamplemap.outputbuffer, resamplemap, A.nodespectra[child])

        _muladd_or_mulreset!(
            _eθ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            conj.(A.phaseshifttoparent[sector, level+1]),
            reset=reset,
        )
        _muladd_or_mulreset!(
            _eϕ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            conj.(A.phaseshifttoparent[sector, level+1]),
            reset=reset,
        )

        reset = false
    end
end

"""
    _disaggregate_to_disaggregationslist!(receivestruct::MLFMMReceive)

Perform all disaggregations according to `receivestruct.disaggregationlist`.
The `receivestruct.disaggregationlist` stores all nodes at each level which shall perform a disaggregation.
"""
function _disaggregate_to_disaggregationslist!(receivestruct::MLFMMReceive)
    receivestruct.verbose && @info "Disaggregate node spectra "


    for level in eachindex(receivestruct.disaggregationlist)
        receivestruct.disaggregationlist[level] == [] && continue
        for sector in 1:8
            # Threads.@threads for sector = 1:8
            conj!(receivestruct.phaseshifttoparent[sector, level+1])
        end

        for node in receivestruct.disaggregationlist[level]
            # Threads.@threads for node in receivestruct.disaggregationlist[level]
            !(receivestruct.nodeisoccupied[node]) && continue
            _disaggregate_children!(receivestruct, node)
        end

        for sector in 1:8
            # Threads.@threads for sector = 1:8
            conj!(receivestruct.phaseshifttoparent[sector, level+1])
        end

    end

    _disaggregate_leafnodes!(receivestruct)

end

"""
    _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive, [y::AbstractVector])

Perform adjoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate!`#

If no vector `y` is given, the content in `A.buffer` is used for adjoint disaggregation. 
Otherwise, `A.buffer` is overwritten by `y` before adjoint disaggregation.
"""
function _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive, y::AbstractVector)
    A.buffer .= y
    _adjoint_disaggregate_to_disaggregationslist!(A)
end
function _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive)
    A.verbose && @info "Adjoint disaggregate node spectra "

    _adjoint_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))

        for node in A.disaggregationlist[level]
            # Threads.@threads for node in A.disaggregationlist[level]
            !(A.nodeisoccupied[node]) && continue
            # _adjoint_disaggregate_children!(A, node)
            _transpose_disaggregate_children!(A, node)
        end

    end
    # _adjoint_disaggregate_leafnodes!(A)

end

"""
    _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive, [y::AbstractVector])

    Perform transposed operation of `_disaggregate!`

If no vector `y` is given, the content in `A.buffer` is used for transposed disaggregation. 
Otherwise, `A.buffer` is overwritten by `y` before transposed disaggregation.
"""
function _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive, y::AbstractVector)
    A.buffer .= y
    _transpose_disaggregate_to_disaggregationslist!(A)
end
function _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive)
    A.verbose && @info "Transpose disaggregate node spectra "

    _transpose_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))
        A.disaggregationlist[level] == [] && continue

        # Threads.@threads for sector = 1:8
        for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end


        # Threads.@threads for node in A.disaggregationlist[level]
        for node in A.disaggregationlist[level]
            !(A.nodeisoccupied[node]) && continue
            _transpose_disaggregate_children!(A, node)
        end

        # Threads.@threads for sector = 1:8
        for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end

    end

end

"""
    _transfer!(A::MLFMMTransmitMap)

Perform all required transfers between the `sourcestruct` and the `receivestruct` of `A`.
"""
function _transfer!(A::MLFMMTransmitMap)
    transferlist = A.transferlist
    transferplan = A.transferplan
    A.verbose && @info "Transfer source → receive"

    sourcestruct = A.sourcestruct
    receivestruct = A.receivestruct
    receivestruct.nodeisfresh .= false


    # Threads.@threads for receivenode in receivestruct.receivetranslationnodes
    for receivenode in A.receivetranslationnodes
        transfers = transferlist[receivenode]
        for transfernode in transfers
            transfer!(
                receivestruct.nodespectra[receivenode],
                sourcestruct.nodefarfields[transfernode],
                transferplan[receivenode][transfernode],
                reset=!receivestruct.nodeisfresh[receivenode],
            )
            receivestruct.nodeisfresh[receivenode] = true
        end
    end
    nothing
end

"""
    _adjoint_transfer!(A::MLFMMTransmitMap)

Perform the adjoint operator (i.e., complex conjugate of transposed operator) of `_transfer!(A::MLFMMTransmitMap)`
"""
function _adjoint_transfer!(A::MLFMMTransmitMap)

    adjoint_transferlist = A.adjoint_transferlist
    transferplan = A.transferplan
    A.verbose && @info "Adjoint transfer source ← receive"

    sourcestruct = A.sourcestruct
    receivestruct = A.receivestruct

    sourcestruct.nodeisfresh .= false

    # Threads.@threads for transfernode in receivestruct.firetranslationnodes
    for transfernode in A.firetranslationnodes
        transfers = adjoint_transferlist[transfernode]
        for receivenode in transfers
            _adjoint_transfer!(
                receivestruct.nodespectra[receivenode],
                sourcestruct.nodefarfields[transfernode],
                transferplan[receivenode][transfernode],
                reset=!sourcestruct.nodeisfresh[transfernode],
            )
            sourcestruct.nodeisfresh[transfernode] = true
        end
    end
    nothing

end

"""
    _transpose_transfer!(A::MLFMMTransmitMap)

Perform the transposed operator to `_transfer!(A::MLFMMTransmitMap)`
"""
function _transpose_transfer!(A::MLFMMTransmitMap)

    adjoint_transferlist = A.adjoint_transferlist
    transferplan = A.transferplan
    A.verbose && @info "transpose transfer source ← receive"

    sourcestruct = A.sourcestruct
    receivestruct = A.receivestruct

    sourcestruct.nodeisfresh .= false

    # Threads.@threads for transfernode in receivestruct.firetranslationnodes
    for transfernode in A.firetranslationnodes
        transfers = adjoint_transferlist[transfernode]
        for receivenode in transfers
            _transpose_transfer!(
                receivestruct.nodespectra[receivenode],
                sourcestruct.nodefarfields[transfernode],
                transferplan[receivenode][transfernode],
                reset=!sourcestruct.nodeisfresh[transfernode],
            )
            sourcestruct.nodeisfresh[transfernode] = true
        end
    end
    nothing

end

"""
    _forward!((A::MLFMMTransmitMap), [x])

Store the result of the matrix vector product in `A.outputbuffer`

If no vector `x` is given, `A.inputbuffer` is used as excitation vector.
Otherwise, `A.inputvector` is overwritten by the content of `x` before the operation is performed
"""
function _forward!(A::MLFMMTransmitMap)
    A.verbose && @info "-------------------------\n   Evaluate forward operator\n-------------------------------"

    # _aggregate_to_minlevel!(sourcestruct; min_aggregationlevel=receivestruct.minsourcetranslationlevel)
    _aggregate_to_aggregationlist!(A.sourcestruct)
    _transfer!(A)
    _disaggregate_to_disaggregationslist!(A.receivestruct)
    A.outputbuffer .= A.receivestruct.buffer

    A.verbose &&
        println("---------------------------------")
    nothing
end
function _forward!(A::MLFMMTransmitMap, x)
    A.inputbuffer .= x
    A.sourcestruct.buffer .= x
    _forward!(A)
end

"""
    _transpose_forward!(A::MLFMMTransmitMap, [y])

Store the result of the transposed matrix vector product in `A.inputbuffer`

If no vector `y` is given, `A.outputbuffer` is used as input vector.
Otherwise, `A.outputbuffer` is overwritten by the content of `y` before the operation is performed
"""
function _transpose_forward!(A::MLFMMTransmitMap)
    A.verbose &&
        @info "---------------------------------\n   Evaluate transpose forward operator\n---------------------------------------"

    _transpose_disaggregate_to_disaggregationslist!(A.receivestruct)
    _transpose_transfer!(A)
    _transpose_aggregate_to_aggregationlist!(A.sourcestruct)
    A.inputbuffer .= A.sourcestruct.buffer

    A.verbose &&
        println("-----------------------------------------")
    nothing
end
function _transpose_forward!(A::MLFMMTransmitMap, y)
    A.outputbuffer .= y
    A.receivestruct.buffer .= y
    _transpose_forward!(A)
end

"""
    _adjoint_forward!(A::MLFMMTransmitMap)

Store the result of the adjoint matrix vector product in `A.inputbuffer`

If no vector `y` is given, `A.outputbuffer` is used as input vector.
Otherwise, `A.outputbuffer` is overwritten by the content of `y` before the operation is performed    
"""
function _adjoint_forward!(A::MLFMMTransmitMap)
    A.verbose &&
        @info "---------------------------------\n   Evaluate adjoint forward operator\n---------------------------------------"

    _adjoint_disaggregate_to_disaggregationslist!(A.receivestruct)
    _adjoint_transfer!(A)
    _adjoint_aggregate_to_aggregationlist!(A.sourcestruct)
    A.inputbuffer .= A.sourcestruct.buffer

    A.verbose &&
        println("-----------------------------------------")
    nothing
end
function _adjoint_forward!(A::MLFMMTransmitMap, y)
    A.outputbuffer .= y
    A.receivestruct.buffer .= y
    _adjoint_forward!(A)
end

function LinearMaps._unsafe_mul!(y, A::MLFMMTransmitMap, x)
    _forward!(A, x)
    y .= A.outputbuffer
end

function LinearMaps._unsafe_mul!(y, A_ad::LinearMaps.AdjointMap{C,MLFMMTransmitMap{AA,F,M,Y,C,X,TP}}, x::AbstractVector) where {AA,F,M,Y,C,X,TP}
    A = A_ad.lmap
    _adjoint_forward!(A, x)
    y .= A.inputbuffer
end

function LinearMaps._unsafe_mul!(y, A_tr::LinearMaps.TransposeMap{C,MLFMMTransmitMap{AA,F,M,Y,C,X,TP}}, x::AbstractVector) where {AA,F,M,Y,C,X,TP}
    A = A_tr.lmap
    _transpose_forward!(A, x)
    y .= A.inputbuffer
end