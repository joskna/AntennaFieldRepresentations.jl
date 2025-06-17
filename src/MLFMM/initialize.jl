
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
    transfertype::P = PlannedTransfer{
        typeof(nodespectra[3].sampling),
        eltype(nodespectra[3].EθEϕ),
        typeof(nodespectra[3].k0),
    },
    verbose::Bool = false,
) where {P<:Type{<:AbstractTransfer}}
    verbose && @info "Initialize   transfers"

    transferlist = [Int[] for _ = 1:numberofnodes(receivetree)]
    adjoint_transferlist = [Int[] for _ = 1:numberofnodes(sourcetree)]
    mintranslationlevel = typemax(Int)
    uniquetransfers = [transfertype[] for _ in levels(receivetree)]

    transferplan = [spzeros(Int, length(sourcetree.nodes)) for _ in eachindex(transferlist)]

    Pℓstorage =
        Vector{typeof(nodespectra[3].wavenumber)}(undef, levelcutoffparameters[1] + 1)


    for receivenode in DepthFirstIterator(receivetree, root(receivetree))
        !(receivenodeisoccupied[receivenode]) && continue

        receivelevel = level(receivetree, receivenode)

        for sourcenode in nodesatlevel(sourcetree, receivelevel)
            !(transmitnodeisoccupied[sourcenode]) && continue

            if _transfercanhappen(
                sourcenode,
                receivenode,
                sourcetree,
                numbufferboxes = numbufferboxes,
            )
                append!(transferlist[receivenode], sourcenode)
                append!(adjoint_transferlist[sourcenode], receivenode)
                mintranslationlevel = minimum([Int(receivelevel), mintranslationlevel])
                boxhalfsize = halfsize(receivetree, receivenode)
                transvector =
                    center(receivetree, receivenode) - center(sourcetree, sourcenode)
                isnew = true
                for (uniqueindex, uniquetransfer) in
                    enumerate(uniquetransfers[receivelevel])
                    uniquetransfervector = gettransfervector(uniquetransfer)
                    if norm(transvector - uniquetransfervector) < 1e-3 * boxhalfsize
                        isnew = false

                        transferplan[receivenode][sourcenode] = uniqueindex
                        break
                    end
                end

                if isnew
                    mintranslationlevel = minimum([mintranslationlevel, receivelevel])

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

                    # transferplan[receivenode][sourcenode] = _initialize_plannedtransfer!(
                    #     Pℓstorage,
                    #     transvector,
                    #     getwavenumber(nodefarfields[sourcenode]),
                    #     nodefarfields[sourcenode].samplingstrategy,
                    #     levelcutoffparameters[receivelevel],
                    #     multiplyweights = true,
                    # )

                    # append!(
                    #     uniquetransfers[receivelevel],
                    #     [transferplan[receivenode][sourcenode]],
                    # )
                    append!(
                        uniquetransfers[receivelevel],
                        [
                            _initialize_plannedtransfer!(
                                Pℓstorage,
                                transvector,
                                getwavenumber(nodefarfields[sourcenode]),
                                nodefarfields[sourcenode].samplingstrategy,
                                levelcutoffparameters[receivelevel],
                                multiplyweights = true,
                            ),
                        ],
                    )
                    transferplan[receivenode][sourcenode] =
                        length(uniquetransfers[receivelevel])

                end
            end

        end

    end

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
function _transfercanhappen(sourcenode, receivenode, tree; numbufferboxes = 1)
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
function _findoccupiednodes(tree; minlevel = 1)
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
    minhalfsize::R;
    verbose = false,
) where {R<:Real}
    verbose && @info "Initialize  MLFMM tree"
    rootcenter, rootsize = getboundingbox(points)

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2 = maximum([0, Integer(ceil(log2(rootsize / minhalfsize)))])
    rootsize = 2^halfsize_powerof2 * minhalfsize

    root_center = SVector{3,Float64}(rootcenter)

    tree = AntennaFieldRepresentations.MLFMMTree(root_center, points, rootsize, minhalfsize)

    message = string(" ", length(levels(tree)), " levels,  ", numberofnodes(tree), " boxes")
    verbose && @info message
    return tree
end
function _initialize_trees(
    sourcepoints::AbstractArray{SVector{3,R}},
    receivepoints::AbstractArray{SVector{3,R}},
    minhalfsize::R;
    verbose = false,
) where {R<:Real}
    verbose && @info "Initialize  MLFMM trees"
    srcrootcenter, srcroothalflength = getboundingbox(sourcepoints)

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2 = maximum([0, Integer(ceil(log2(srcroothalflength / minhalfsize)))])
    srcroothalflength = 2^halfsize_powerof2 * minhalfsize

    # ensure that center is at one corner of sourcebox
    sourceroot_center = SVector{3,Float64}(
        srcrootcenter - [srcroothalflength, srcroothalflength, srcroothalflength],
    )

    recroot_center, recroothalflength = getboundingbox(receivepoints)

    roothalflength = maximum([
        srcroothalflength,
        recroothalflength + maximum(abs.(sourceroot_center - recroot_center)),
    ])

    #ensure that halfsize is 2ᴺ ⋅ srcroothalflength
    halfsize_powerof2 =
        maximum([0, Integer(ceil(log2(roothalflength / srcroothalflength)))])
    roothalflength = 2^halfsize_powerof2 * srcroothalflength

    points = [sourcepoints; receivepoints]
    tree = AntennaFieldRepresentations.MLFMMTree(
        sourceroot_center,
        points,
        roothalflength,
        minhalfsize,
    )

    message = string(" ", length(levels(tree)), " levels,  ", numberofnodes(tree), " boxes")
    verbose && @info message
    return tree
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
    minlevel::Integer = 1,
    verbose::Bool = false,
) where {P<:PropagationType,S<:SphereSamplingStrategy,C}
    verbose && println()
    verbose && @info string("Allocate node patterns of type ", W)
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
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
    verbose::Bool = false,
) where {T<:Real,S,Y<:SphereSamplingStrategy}
    verbose && @info "Compute basis patterns"
    L = cutoffparameter
    samplingstrategy = _standardsampling(samplingtype, L)
    return _initializebasisfunctionfarfields(
        tree,
        basisfunctions,
        samplingstrategy,
        k0;
        verbose = verbose,
    )

end
function _initializebasisfunctionfarfields(
    tree::MLFMMTree,
    basisfunctions::S,
    samplingstrategy::Y,
    k0::T;
    verbose::Bool = false,
) where {T<:Real,S,Y<:SphereSamplingStrategy}
    verbose && @info "Compute basis far fields"
    message = string(length(basisfunctions), " basis functions")
    verbose && @info message

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
    samplingstrategy::Y;
    verbose::Bool = false,
) where {IFS<:IrregularFieldSampling,Y<:SphereSamplingStrategy}
    verbose && @info "Compute probe patterns"
    message = string(length(fieldsampling), " sampling positions")
    verbose && @info message
    θvec, ϕvec = samples(samplingstrategy)
    θvec .= pi .- θvec
    ϕvec .= mod2pi.(ϕvec .+ pi)

    nθ, nϕ = length(θvec), length(ϕvec)

    C = eltype(fieldsampling)
    positions = fieldsampling.positions


    k0 = getwavenumber(fieldsampling.probes[1].aut_field)


    basisfunctionweightingpatterns =
        Array{PlaneWaveExpansion{Radiated,Y,C}}(undef, size(fieldsampling.probeIDs))


    for k in eachindex(basisfunctionweightingpatterns)
        χ, θ, ϕ = fieldsampling.eulerangles[k]

        probefield =
            rotate(fieldsampling.probes[fieldsampling.probeIDs[k]].aut_field, χ, θ, ϕ)
        basisfunctionweightingpatterns[k] = PlaneWaveExpansion(
            Radiated(),
            samplingstrategy,
            zeros(C, nθ, nϕ),
            zeros(C, nθ, nϕ),
            k0,
        )

        for (θind, θangle) in enumerate(θvec)
            for (ϕind, ϕangle) in enumerate(ϕvec)
                _eθ(basisfunctionweightingpatterns[k])[θind, ϕind],
                _eϕ(basisfunctionweightingpatterns[k])[θind, ϕind] =
                    farfield(probefield, (θangle, ϕangle))

                _eϕ(basisfunctionweightingpatterns[k])[θind, ϕind] *= -1
            end
        end

    end


    phaseshiftmatrix = Array{C}(undef, nθ, nϕ)
    for leafnode::Int in leafs(tree)
        for functionindex::Int in tree(leafnode).data.values
            phaseshiftmatrix .= _phaseshiftmatrix!(
                phaseshiftmatrix,
                positions[functionindex] - center(tree, leafnode),
                k0,
                samplingstrategy,
            )
            _eθ(basisfunctionweightingpatterns[functionindex]) .=
                _eθ(basisfunctionweightingpatterns[functionindex]) .* phaseshiftmatrix
            _eϕ(basisfunctionweightingpatterns[functionindex]) .=
                _eϕ(basisfunctionweightingpatterns[functionindex]) .* phaseshiftmatrix
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
    minlevel::Int = 0,
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
    verbose::Bool = false,
) where {Y}
    verbose && @info "Assemble resample maps"
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
    tree::MLFMMTree,
    cutoffparameters::Vector{I},
    k0::T;
    minlevel::Int = 0,
    samplingtype::Type{Y} = GaussLegendreθRegularϕSampling,
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
