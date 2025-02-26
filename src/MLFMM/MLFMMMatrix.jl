include("beastglue.jl")
include("MLFMMTree.jl")

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
    leafnodeindices::Vector{Int}
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
function setwavenumber!(p::MLFMMSource{M,Y,C,X}, val) where {M,Y,C,X}
    p = MLFMMSource{M,Y,C,X}(
        p.tree,
        p.expectedaccuracy,
        val,
        p.nodefarfields,
        p.basisfunctionfarfields,
        p.levelcutoffparameters,
        p.levelresamplemaps,
        p.phaseshifttoparent,
        p.nodeisfresh,
        p.leafnodeindices,
        p.buffer,
        p.outputbuffer,
        p.verbose,
    )
    return p
end


struct MLFMMTransmitMap{
    A<:AntennaFieldRepresentation,
    F<:FieldSampling,
    M<:ResampleMap,
    Y<:SphereSamplingStrategy,
    C<:Complex,
} <: TransmitMap{A,F,C}
    inputbuffer::Vector{C}
    outputbuffer::Vector{C}
    tree::MLFMMTree
    expectedaccuracy::Real
    wavenumber::Real
    nodefarfields::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    nodespectra::Vector{PlaneWaveExpansion{Incident,Y,C}}
    basisfunctionfarfields::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    testfunctionfarfieldds::Vector{PlaneWaveExpansion{Radiated,Y,C}}
    levelcutoffparameters::Vector{Int}
    levelresamplemaps::Vector{M}
    phaseshifttoparent::Array{Matrix{C},2}
    transferlist::Vector{Vector{Int}}
    adjoint_transferlist::Vector{Vector{Int}}
    transferplan::Vector{Vector{PlannedTransfer{C}}}
    transmitnodeisfresh::Vector{Bool}
    receivenodeisfresh::Vector{Bool}
    firetranslationnodes::Vector{Int}
    aggregationlist::Vector{Vector{Int}}
    disaggregationlist::Vector{Vector{Int}}
    transmitleafnodeindices::Vector{Int}
    receiveleafnodeindices::Vector{Int}
    verbose::Bool
    tmpmatrix::Matrix{C}
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

    nodefarfields = _allocatenodepattern(
        PlaneWaveExpansion{Radiated,samplingtype,Complex{T}},
        tree,
        levelcutoffparameters,
        wavenumber,
    )
    outputbuffer = nodefarfields[1].buffer

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
        T,
        tree,
        levelcutoffparameters,
        wavenumber;
        samplingtype = samplingtype,
    )
    L = levelcutoffparameters[1]

    sampling = _standardsampling(samplingtype, L)
    θs, ϕs = samples(sampling)
    nθ, nϕ = length(θs), length(ϕs)
    R = center(tree, 1)
    globalphaseshift = Matrix{Complex{T}}(undef, nθ, nϕ)
    _phaseshiftmatrix!(globalphaseshift, -R, wavenumber, sampling)

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
        leafnodeindices,
        inputbuffer,
        outputbuffer,
        verbose,
    )

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
function _allocatenodepattern(
    W::Type{PlaneWaveExpansion{P,S,C}},
    tree::MLFMMTree,
    cutoffparameters::Vector{<:Integer},
    wavenumber;
    minlevel::Integer = 1,
) where {P<:PropagationType,S<:SphereSamplingStrategy,C}
    nodepattern = Vector{W}(undef, length(tree.nodes))

    levels = AntennaFieldRepresentations.levels(tree)
    for level = minlevel:length(levels)
        L = cutoffparameters[level]
        samplingstrategy = _standardsampling(S, L)
        θs, ϕs = samples(samplingstrategy)
        for node in AntennaFieldRepresentations.nodesatlevel(tree, level)
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
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
    verbose::Bool = false,
) where {T<:Real,S,Y<:SphereSamplingStrategy}
    L = cutoffparameter
    samplingstrategy = _standardsampling(samplingtype, L)
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
    basisfunctions::SurfaceCurrentDensity,
    sampling::Y,
) where {Y<:SphereSamplingStrategy}
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
        point(cosp[kk] * sint[k], sinp[kk] * sint[k], cost[k]) for k in eachindex(θvec)
        for kk in eachindex(ϕvec)
    ]

    farfieldmatrices = individualcartesianfarfields(
        basisfunctions,
        pts,
    )::Matrix{SVector{3,Complex{eltype(θvec)}}}

    basisfunctionfarfields = Vector{PlaneWaveExpansion{Radiated,Y,Complex{T}}}(
        undef,
        numfunctions(basisfunctions),
    )
    ff = Matrix{SVector{3,Complex{T}}}(undef, nθ, nϕ)
    for functionindex in eachindex(basisfunctionfarfields)
        ff .= reshape(view(farfieldmatrices, :, functionindex), nθ, nϕ)

        basisfunctionfarfields[functionindex] = PlaneWaveExpansion(
            Radiated(),
            sampling,
            Matrix{Complex{T}}(undef, nθ, nϕ),
            Matrix{Complex{T}}(undef, nθ, nϕ),
            getwavenumber(basisfunctions),
        )
        Eθ = _eθ(basisfunctionfarfields[functionindex])
        Eϕ = _eϕ(basisfunctionfarfields[functionindex])
        for kθ in eachindex(θvec), kϕ in eachindex(ϕvec)


            Eθ[kθ, kϕ] = Complex{R}(udot(eθ[kθ, kϕ], ff[kθ, kϕ]))
            Eϕ[kθ, kϕ] = Complex{R}(udot(eϕ[kϕ], ff[kθ, kϕ]))
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
    _initializelevelresamplemaps(R::Type{<:Real}, tree::MLFMMTree, orderθ::Integer, orderϕ::Integer,cutoffparameters::Vector{<:Integer}; minlevel::Int=0)

Return list of interpolators to aggregate from each level to the respective parent level.
"""
function _initializelevelresamplemaps(
    T::Type{<:Real},
    tree::MLFMMTree,
    orderθ::Integer,
    orderϕ::Integer,
    cutoffparameters::Vector{<:Integer};
    minlevel::Int = 0,
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
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
            orderθ = orderθ,
            orderϕ = orderϕ,
        )
    end
    return levelresamplemaps
end
#TODO: only initialize required phaseshifts and share between source and receive
function _initializephaseshifttoparent(
    T::Type{<:Real},
    tree::MLFMMTree,
    cutoffparameters::Vector{<:Integer},
    k0::Real;
    minlevel::Int = 0,
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
) where {Y}
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

    # for leafnode in A.leafnodeindices
    Threads.@threads for leafnode in A.leafnodeindices
        reset = true
        for functionindex::Int in tree(leafnode).data.values::Vector{Int}
            A.nodefarfields[leafnode] .= _muladd_or_mulreset!(
                A.nodefarfields[leafnode],
                A.basisfunctionfarfields[functionindex],
                A.buffer[functionindex],
                reset = reset,
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
    Threads.@threads for leafnode in A.leafnodeindices
        pws = A.nodefarfields[leafnode]
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff = A.basisfunctionfarfields[basisfunctionindex]
            A.xvector[basisfunctionindex] = dot(ff, pws)
        end
    end
end
"""
    _transpose_aggregate_leafnodes!(A::MLFMMSource)

Perform the transpose operation to `_aggregate_leafnodes!`
"""
function _transpose_aggregate_leafnodes!(A::MLFMMSource)
    tree = A.tree
    Threads.@threads for leafnode in A.leafnodeindices
        pws = A.nodefarfields[leafnode]
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff = A.basisfunctionfarfields[basisfunctionindex]
            A.xvector[basisfunctionindex] = udot(ff, pws)
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

        sector = tree.nodes[child].data.sector + 1

        resamplemap.outputbuffer .=
            mul!(resamplemap.outputbuffer, resamplemap, A.nodefarfields[child])
        _eθ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(
            _eθ(A.nodefarfields[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            A.phaseshifttoparent[sector, level+1],
            reset = reset,
        )
        _eϕ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(
            _eϕ(A.nodefarfields[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            A.phaseshifttoparent[sector, level+1],
            reset = reset,
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

        sector = tree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        A.nodefarfields[child] .=
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

        sector = tree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(A.nodefarfields[parentnode]) .* conj.(A.phaseshifttoparent[sector, level+1])
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(A.nodefarfields[parentnode]) .* conj.(A.phaseshifttoparent[sector, level+1])
        A.nodefarfields[child] .=
            mul!(A.nodefarfields[child], adjoint_resamplemat, resamplemap.outputbuffer)


        # resamplemap.outputbuffer .= mul!(resamplemap.outputbuffer, resamplemap, A.nodefarfields[child])
        # _eθ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eθ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 1),  A.phaseshifttoparent[sector, level+1], reset = reset)
        # _eϕ(A.nodefarfields[parentnode]) .= _muladd_or_mulreset!(_eϕ(A.nodefarfields[parentnode]), view(resamplemap.outputbuffermat, :, :, 2),  A.phaseshifttoparent[sector, level+1], reset = reset)

        # reset = false
    end
end

"""
    _aggregate_to_minlevel!(A::MLFMMSource, [x::AbstractVector]; min_aggregationlevel::Integer=0)

Aggregate `A` up to min_aggregationlevel. 
"""
function _aggregate_to_minlevel!(A::MLFMMSource, x; min_aggregationlevel::Integer = 0)
    A.buffer .= x
    _aggregate_to_minlevel!(A, min_aggregationlevel = min_aggregationlevel)
end
function _aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer = 0)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    _aggregate_leafnodes!(A)

    levels = AntennaFieldRepresentations.levels(tree)
    for level in reverse(maximum([min_aggregationlevel, 1]):length(levels))
        for parentnode::Int in nodesatlevel(tree, level)
            isleaf(tree, parentnode) && continue
            _aggregate_children!(A, parentnode)
        end
    end

end
"""
    _adjoint_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)


Perform adjoint operation to `aggregate_to_minlevel!`
"""
function _adjoint_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer = 0)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    levels = AntennaFieldRepresentations.levels(tree)
    for level = maximum([min_aggregationlevel, 1]):length(levels)
        # for parentnode::Int in nodesatlevel(tree, level)
        for parentnode::Int in nodesatlevel(tree, level)
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
    min_aggregationlevel::Integer = 0,
)
    A.verbose && @info "Aggregate node far fields"
    tree = A.tree

    levels = AntennaFieldRepresentations.levels(tree)
    for level = maximum([min_aggregationlevel, 1]):length(levels)
        # for parentnode::Int in nodesatlevel(tree, level)
        for parentnode::Int in nodesatlevel(tree, level)
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
    _eθ(A.nodefarfields[1]) .= _eθ(A.nodefarfields[1]) .* A.globalphaseshift
    _eϕ(A.nodefarfields[1]) .= _eϕ(A.nodefarfields[1]) .* A.globalphaseshift
end
"""
    _adjoint_aggregate_to_farfield!(A::MLFMMSource; min_aggregationlevel::Integer=0)

Perform adjoint operation to `_aggregate_to_farfield!`
"""
function _adjoint_aggregate_to_farfield!(A::MLFMMSource)
    _eθ(A.nodefarfields[1]) .= _eθ(A.nodefarfields[1]) .* conj.(A.globalphaseshift)
    _eϕ(A.nodefarfields[1]) .= _eϕ(A.nodefarfields[1]) .* conj.(A.globalphaseshift)
    _adjoint_aggregate_to_minlevel!(A)
end
"""
    _transpose_aggregate_to_farfield!(A::MLFMMSource; min_aggregationlevel::Integer=0)

Perform transpose operation to `_aggregate_to_farfield!`
"""
function _transpose_aggregate_to_farfield!(A::MLFMMSource)
    _eθ(A.nodefarfields[1]) .= _eθ(A.nodefarfields[1]) .* A.globalphaseshift
    _eϕ(A.nodefarfields[1]) .= _eϕ(A.nodefarfields[1]) .* A.globalphaseshift
    _transpose_aggregate_to_minlevel!(A)
end


# struct MLFMMReceive{R<:ResampleMap,Y<:SphereSamplingStrategy,R<:Real}
#     tree::MLFMMtree
#     # expectedaccuracy::R
#     # wavenumber::R
#     nodespectra::Vector{PlaneWaveExpansion{Incident,Y,Complex{R}}}
#     basisfunctionpatterns::Vector{PlaneWaveExpansion{Radiated,Y,Complex{R}}}
#     levelcutoffparameters::Vector{Int}
#     levelinterpolators::Vector{A}
#     phaseshifttoparent::Array{Matrix{Complex{R}},2}
#     transferlist::Vector{Vector{Int}}
#     adjoint_transferlist::Vector{Vector{Int}}
#     transferplan::Vector{Vector{PlannedTransfer{Complex{R}}}}
#     nodeisfresh::Vector{Bool}
#     leafnodeindices::Vector{Int}
#     minreceivelevel::Int
#     minsourcetranslationlevel::Int
#     receivetranslationnodes::Vector{Int}
#     firetranslationnodes::Vector{Int}
#     aggregationlist::Vector{Vector{Int}}
#     disaggregationlist::Vector{Vector{Int}}
#     bvector::Vector{Complex{R}}
#     verbose::Bool
#     tmpmatrix::Matrix{Complex{R}}
# end
