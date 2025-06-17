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
    expectedaccuracy = T(1e-3),
    verbose = false,
    minhalfsize = π / (2 * wavenumber),
    orderθ = 8,
    orderϕ = 8,
    samplingtype::Type{S} = GaussLegendreθRegularϕSampling,
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
        verbose = verbose,
        samplingtype = samplingtype,
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
        samplingtype = samplingtype,
    )
    phaseshifttoparent = _initializephaseshifttoparent(
        tree,
        levelcutoffparameters,
        T(wavenumber);
        samplingtype = samplingtype,
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
    transferplan::Vector{SparseVector{Int,Int}}
    uniquetransfers::Vector{Vector{TP}}
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
    expectedaccuracy = T(1e-3),
    verbose = false,
    minhalfsize = π / (2 * wavenumber),
    orderθ = 8,
    orderϕ = 8,
    samplingtype::Type{S} = GaussLegendreθRegularϕSampling,
    num_bufferboxes::Integer = 1,
    transfertype = PlannedTransfer{samplingtype,Complex{T},T},
) where {T<:Real,S<:SphereSamplingStrategy}
    C = Complex{T}

    sourcepoints = _getpoints(basisfunctions)
    receivepoints = _getpoints(fieldsampling)
    points = [sourcepoints; receivepoints]
    # sourcetree = _initialize_tree(points, minhalfsize, verbose=verbose)
    sourcetree =
        _initialize_trees(sourcepoints, receivepoints, minhalfsize, verbose = verbose)

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
    message = string(sum(transmitnodeisoccupied), " boxes contain sources")
    verbose && @info message

    receivenodeisoccupied = _findoccupiednodes(receivetree)
    message = string(sum(receivenodeisoccupied), " boxes contain probes")
    verbose && @info message

    sourcerootnode = _findrootnode(sourcetree, transmitnodeisoccupied)
    minsourcelevel = level(sourcetree, sourcerootnode)
    message = string(
        "Source root node ",
        sourcerootnode,
        " at level ",
        level(sourcetree, sourcerootnode),
        " contains all sources",
    )
    verbose && @info message

    receiverootnode = _findrootnode(receivetree, receivenodeisoccupied)
    minreceivelevel = level(receivetree, receiverootnode)
    message = string(
        "Receive root node ",
        receiverootnode,
        " at level ",
        level(receivetree, receiverootnode),
        " contains all probes",
    )
    verbose && @info message



    transmitnodeisfresh = [false for k = 1:length(sourcetree.nodes)]
    receivenodeisfresh = [false for k = 1:length(receivetree.nodes)]


    nodefarfields = _allocatenodepattern(
        PlaneWaveExpansion{Radiated,samplingtype,Complex{T}},
        sourcetree,
        levelcutoffparameters,
        transmitnodeisoccupied,
        wavenumber,
        minlevel = minimum([3, minsourcelevel]),
        verbose = verbose,
    )

    nodespectra = _allocatenodepattern(
        PlaneWaveExpansion{Incident,samplingtype,Complex{T}},
        receivetree,
        levelcutoffparameters,
        receivenodeisoccupied,
        wavenumber,
        minlevel = 3,
        verbose = verbose,
    )


    L = levelcutoffparameters[end]
    samplingstrategy = _standardsampling(S, L)

    basisfunctionfarfields = _initializebasisfunctionfarfields(
        sourcetree,
        basisfunctions,
        samplingstrategy,
        wavenumber;
        verbose = verbose,
    )
    inputbuffer = Vector{Complex{T}}(undef, length(basisfunctions))
    inputbuffer .= basisfunctions

    # message = string(Base.format_bytes(memory_basispatterns), " in memory")
    # verbose && @info message

    outputbuffer = deepcopy(asvector(fieldsampling))
    # verbose && println()


    testfunctionfarfields = _initializebasisfunctionweightingpatterns(
        receivetree,
        fieldsampling,
        samplingstrategy,
        verbose = verbose,
    )


    levelresamplemaps = _initializelevelresamplemaps(
        T,
        sourcetree,
        orderθ,
        orderϕ,
        levelcutoffparameters;
        samplingtype = samplingtype,
        verbose = verbose,
    )
    phaseshifttoparent = _initializephaseshifttoparent(
        sourcetree,
        levelcutoffparameters,
        T(wavenumber);
        samplingtype = samplingtype,
    )



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
        transfertype = transfertype,
        verbose = verbose,
    )


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
            outputbuffer,
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
        uniquetransfers,
        firetranslationnodes,
        mintranslationlevel,
        receivetranslationnodes,
        verbose,
    )

    verbose && _show_memory(transmap)

    verbose && println("------------------------------")

    return transmap

end

function _show_memory(A::MLFMMTransmitMap)
    memory_tree =
        Base.summarysize(A.sourcestruct.tree) + Base.summarysize(A.receivestruct.tree)
    memory_nodefarfields = Base.summarysize(A.sourcestruct.nodefarfields)
    memory_basispatterns = Base.summarysize(A.sourcestruct.basisfunctionfarfields)
    memory_testfunctionfarfields = Base.summarysize(A.receivestruct.testfunctionfarfields)
    memory_nodesprectra = Base.summarysize(A.receivestruct.nodespectra)
    memory_resamplemaps = Base.summarysize(A.sourcestruct.levelresamplemaps)
    memory_phaseshifttoparent = Base.summarysize(A.sourcestruct.phaseshifttoparent)

    memory_transfers = Base.summarysize(A.uniquetransfers)
    memorylists =
        Base.summarysize(A.sourcestruct.aggregationlist) +
        Base.summarysize(A.receivestruct.disaggregationlist) +
        Base.summarysize(A.transferlist) +
        Base.summarysize(A.adjoint_transferlist)

    memory_globalphaseshift = Base.summarysize(A.sourcestruct.globalphaseshift)
    memorylists += Base.summarysize(A.transferplan)

    memory_total =
        memory_globalphaseshift +
        memory_tree +
        memory_nodesprectra +
        memory_basispatterns +
        memory_testfunctionfarfields +
        memory_nodefarfields +
        memory_resamplemaps +
        memory_phaseshifttoparent +
        memory_transfers +
        memorylists
    exact_memory = Base.summarysize(A)

    println()
    println("-----   Used Memory  ---------")
    message = string(
        "    Exact  total   memory: ",
        Base.format_bytes(exact_memory),
        " ( ",
        exact_memory,
        " bytes )",
    )
    @info message
    message = string(
        "Estimated  total   memory: ",
        Base.format_bytes(memory_total),
        " ( ",
        memory_total,
        " bytes, ",
        round(memory_total / exact_memory * 100, digits = 2),
        " % )",
    )
    @info message
    message = string(
        "Memory  for  tree  struct: ",
        Base.format_bytes(memory_tree),
        " ( ",
        round(memory_tree / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for node farfields: ",
        Base.format_bytes(memory_nodefarfields),
        " ( ",
        round(memory_nodefarfields / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory  for node  spectra: ",
        Base.format_bytes(memory_nodesprectra),
        " ( ",
        round(memory_nodesprectra / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for basis patterns: ",
        Base.format_bytes(memory_basispatterns),
        " ( ",
        round(memory_basispatterns / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for test functions: ",
        Base.format_bytes(memory_testfunctionfarfields),
        " ( ",
        round(memory_testfunctionfarfields / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for  resample maps: ",
        Base.format_bytes(memory_resamplemaps),
        " ( ",
        round(memory_resamplemaps / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for phase to parent: ",
        Base.format_bytes(memory_phaseshifttoparent),
        " ( ",
        round(memory_phaseshifttoparent / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory for translation ops: ",
        Base.format_bytes(memory_transfers),
        " ( ",
        round(memory_transfers / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory  for  global  phase: ",
        Base.format_bytes(memory_globalphaseshift),
        " ( ",
        round(memory_globalphaseshift / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
    message = string(
        "Memory  for  control  flow: ",
        Base.format_bytes(memorylists),
        " ( ",
        round(memorylists / memory_total * 100, digits = 2),
        " % of total )",
    )
    @info message
end

include("initialize.jl")
include("aggregate.jl")
include("disaggregate.jl")

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
    # for nodespectrum in receivestruct.nodespectra[receivestruct.nodeisoccupied]
    #     nodespectrum.EθEϕ .*= 0
    # end

    # Threads.@threads for receivenode in receivestruct.receivetranslationnodes
    # for receivenode in A.receivetranslationnodes
    for receivenode in A.receivetranslationnodes
        transfers = transferlist[receivenode]
        for transfernode in transfers
            transfer!(
                receivestruct.nodespectra[receivenode],
                sourcestruct.nodefarfields[transfernode],
                A.uniquetransfers[level(A.sourcestruct.tree, transfernode)][transferplan[receivenode][transfernode]],
                reset = !receivestruct.nodeisfresh[receivenode],
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
                A.uniquetransfers[level(A.sourcestruct.tree, transfernode)][transferplan[receivenode][transfernode]],
                reset = !sourcestruct.nodeisfresh[transfernode],
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
                A.uniquetransfers[level(A.sourcestruct.tree, transfernode)][transferplan[receivenode][transfernode]],
                reset = !sourcestruct.nodeisfresh[transfernode],
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
    A.verbose &&
        @info "-------------------------\n   Evaluate forward operator\n-------------------------------"

    # _aggregate_to_minlevel!(sourcestruct; min_aggregationlevel=receivestruct.minsourcetranslationlevel)
    _aggregate_to_aggregationlist!(A.sourcestruct)
    _transfer!(A)
    _disaggregate_to_disaggregationslist!(A.receivestruct)
    A.outputbuffer .= A.receivestruct.buffer

    A.verbose && println("---------------------------------")
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

    A.verbose && println("-----------------------------------------")
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

    A.verbose && println("-----------------------------------------")
    nothing
end
function _adjoint_forward!(A::MLFMMTransmitMap, y)
    A.outputbuffer .= y
    A.receivestruct.buffer .= y
    _adjoint_forward!(A)
end

function LinearMaps._unsafe_mul!(y, A::MLFMMTransmitMap, x::AbstractVector)
    _forward!(A, x)
    y .= A.outputbuffer
end

function LinearMaps._unsafe_mul!(
    y,
    A_ad::LinearMaps.AdjointMap{C,MLFMMTransmitMap{AA,F,M,Y,C,X,TP}},
    x::AbstractVector,
) where {AA,F,M,Y,C,X,TP}
    A = A_ad.lmap
    _adjoint_forward!(A, x)
    y .= A.inputbuffer
end

function LinearMaps._unsafe_mul!(
    y,
    A_tr::LinearMaps.TransposeMap{C,MLFMMTransmitMap{AA,F,M,Y,C,X,TP}},
    x::AbstractVector,
) where {AA,F,M,Y,C,X,TP}
    A = A_tr.lmap
    _transpose_forward!(A, x)
    y .= A.inputbuffer
end

readablesize(x) = Base.format_bytes(Base.summarysize(x))
