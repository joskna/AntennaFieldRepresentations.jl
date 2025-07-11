# TODO: docstrings
# TODO: Keep original representation in inputbuffer for later reference
struct MLFMMSource{
    A<:AntennaFieldRepresentation,
    M<:ResampleMap,
    Y<:SphereSamplingStrategy,
    C,
    X<:MLFMMTree,
} <: AntennaFieldRepresentation{Radiated,C}
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
    buffer::A
    outputbuffer::Vector{C}
    verbose::Bool
end
function asvector(p::MLFMMSource)
    return asvector(p.buffer)
end
Base.size(p::MLFMMSource) = size(asvector(p))
Base.getindex(p::MLFMMSource, i) = Base.getindex(asvector(p), i)
Base.setindex!(p::MLFMMSource, i, v) = Base.setindex!(asvector(p), i, v)
function Base.copy(p::MLFMMSource{A,M,Y,C,X}) where {A,M,Y,C,X}
    return deepcopy(p)
end
function rotate(p::MLFMMSource, χ, θ, ϕ; orderθ::Integer=12, orderϕ::Integer=12,)
    rotated_aut = rotate(p, χ, θ, ϕ; orderθ=orderθ, orderϕ=orderϕ,)

    MLFMMSource(
        rotated_aut,
        p.wavenumber;
        expectedaccuracy=p.expectedaccuracy,
        verbose=p.verbose,
        minhalfsize=p.tree(leafs(p.tree)[1]).minhalfsize,
        orderθ=p.orderθ,
        orderϕ=p.orderϕ,
        samplingtype=p.samplingtype,
    )
end

# function Base.similar(p::MLFMMSource)
#     return deepcopy(p)
# end

function changerepresentation(
    ::Type{MLFMMSource},
    aut_field::DipoleArray; #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions) 
    wavenumber=getwavenumber(aut_field),
    expectedaccuracy=1e-3,
    verbose=false,
    minhalfsize=π / (2 * wavenumber),
    orderθ=8,
    orderϕ=8,
    samplingtype=GaussLegendreθRegularϕSampling,
)
    return MLFMMSource(
        aut_field,
        wavenumber;
        expectedaccuracy=expectedaccuracy,
        verbose=verbose,
        minhalfsize=minhalfsize,
        orderθ=orderθ,
        orderϕ=orderϕ,
        samplingtype=samplingtype,
    )
end

function changerepresentation(
    ::Type{MLFMMSource},
    aut_field::SurfaceCurrentDensity; #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions) 
    wavenumber=getwavenumber(aut_field),
    expectedaccuracy=1e-3,
    verbose=false,
    minhalfsize=π / (2 * wavenumber),
    orderθ=8,
    orderϕ=8,
    samplingtype=GaussLegendreθRegularϕSampling,
)
    return MLFMMSource(
        aut_field,
        wavenumber;
        expectedaccuracy=expectedaccuracy,
        verbose=verbose,
        minhalfsize=minhalfsize,
        orderθ=orderθ,
        orderϕ=orderϕ,
        samplingtype=samplingtype,
    )
end

function MLFMMSource(
    basisfunctions; #Can be ::SurfaceCurrentDensity, ::DipoleArray, NamedTuple(:points,:sourcefunctions) 
    wavenumber=getwavenumber(basisfunctions),
    expectedaccuracy=1e-3,
    verbose=false,
    minhalfsize=π / (2 * wavenumber),
    orderθ=8,
    orderϕ=8,
    samplingtype=GaussLegendreθRegularϕSampling,
)
    return MLFMMSource(
        basisfunctions,
        wavenumber;
        expectedaccuracy=expectedaccuracy,
        verbose=verbose,
        minhalfsize=minhalfsize,
        orderθ=orderθ,
        orderϕ=orderϕ,
        samplingtype=samplingtype,
    )
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
    # inputbuffer = Vector{Complex{T}}(undef, length(basisfunctions))
    inputbuffer = similar(basisfunctions)
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

    return MLFMMSource{
        typeof(basisfunctions),
        eltype(levelresamplemaps),
        samplingtype,
        Complex{T},
        typeof(tree),
    }(
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
