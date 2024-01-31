#TODO: Implement EFIE and test functions at source funtion location
#TODO: Combined Constructor for Source and Receive
struct MLFMMSource{T<:MLFMMTrees.AbstractMLFMMTree, R<:Real, A<:AbstractInterpolationStrategy}
    tree::T
    expectedaccuracy::R
    wavenumber::R
    nodefarfields::Vector{FarfieldPattern{Complex{R}}}
    basisfunctionfarfields::Vector{FarfieldPattern{Complex{R}}}
    levelcutoffparameters::Vector{Int}
    levelinterpolators::Vector{A}
    phaseshifttoparent::Array{Matrix{Complex{R}},2}
    nodeisfresh::Vector{Bool}
    leafnodeindices::Vector{Int}
    xvector::Vector{Complex{R}}
    verbose::Bool
    tmpmatrix::Matrix{Complex{R}}
end
#TODO: All current types
#TODO: get rid of tree
#TODO: switch between single core and multicore
#TODO: enable and test Float32
function MLFMMSource(basisfunctions, #Can be ::AbstractSurfaceCurrentDensity, ::AbstractArray{D} where{D<:AbstractDipole}, NamedTuple(:points,:sourcefunctions) 
    wavenumber::R;
    expectedaccuracy=R(1e-3),
    verbose=false,
    minhalfsize=π/(2*wavenumber),
    θstencilsize=8,
    ϕstencilsize=8
    ) where {R<:Real}

verbose && @info "----------------------\n   Assemble MLFMM source \n----------------------------"

verbose && @info "Initialize source tree"
points = _getpoints(basisfunctions)
tree = _initialize_tree(points, minhalfsize)

verbose && @info "Allocate node patterns"
levelcutoffparameters = _initializelevelcutoffparameters(tree, expectedaccuracy, wavenumber)
nodefarfields = _allocatenodepattern(FarfieldPattern{Complex{R}}, tree, levelcutoffparameters)
nodeisfresh=[false for _ in 1:length(MLFMMTrees.tree(tree).nodes)]
tmpmatrix=similar(nodefarfields[end].Eθ)

verbose && @info "Assemble leaf patterns"
leafnodeindices = MLFMMTrees.leafs(tree)
basisfunctionfarfields = _initializebasisfunctionfarfields(tree, basisfunctions, levelcutoffparameters, wavenumber, verbose)
xvector=Vector{Complex{R}}(undef, length(basisfunctionfarfields))

verbose && @info "Assemble interpolators"
levelinterpolators = _initializelevelinterpolators(R, tree, θstencilsize, ϕstencilsize, levelcutoffparameters)
phaseshifttoparent = _initializephaseshifttoparent(R, tree,levelcutoffparameters, wavenumber)

verbose && println("------------------------------")
return MLFMMSource(tree, expectedaccuracy, wavenumber, nodefarfields, basisfunctionfarfields, levelcutoffparameters, levelinterpolators, phaseshifttoparent, nodeisfresh, leafnodeindices, xvector, verbose, tmpmatrix)
end

struct MLFMMReceive{T<:MLFMMTrees.AbstractMLFMMTree, R<:Real, A<:AbstractInterpolationStrategy} 
    tree::T
    # expectedaccuracy::R
    # wavenumber::R
    nodespectra::Vector{PlaneWaveSpectrum{Complex{R}}}
    basisfunctionpatterns::Vector{FarfieldPattern{Complex{R}}}
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
#TODO: Share more memory
#TODO: Only acquire node farfields for two levels (previous levels can be overwritten in aggregation and disaggregation) 
function MLFMMReceive(testfunctions, #Can be ::AbstractSurfaceCurrentDensity, ::AbstractArray{D} where{D<:AbstractDipole}, NamedTuple(:points,:sourcefunctions)
    sourcestruct::MLFMMSource;
    verbose=false,
    num_bufferboxes::Integer=1,
    θstencilsize=8,
    ϕstencilsize=8,
    R::DataType = typeof(sourcestruct.wavenumber)
   )
   
wavenumber= sourcestruct.wavenumber
expectedaccuracy=sourcestruct.expectedaccuracy
minhalfsize = MLFMMTrees.tree(sourcestruct.tree).nodes[end].data.halfsize
sourcetree=sourcestruct.tree

verbose && @info "----------------------\n   Assemble MLFMM receive \n----------------------------"
verbose && @info "Initialize  probe tree"
points = _getpoints(testfunctions)
receivetree = _initialize_tree(points, minhalfsize)

verbose && @info "Initialize   transfers"
transferlist, adjoint_transferlist, minreceivelevel, minsourcetranslationlevel, receivetranslationnodes, firetranslationnodes = _initialize_transferlist(sourcetree, receivetree, num_bufferboxes=num_bufferboxes)
transferplan = _initializetransferplan(transferlist, sourcestruct, receivetree; transfertype=PlannedTransfer{Complex{R}})
nodeisfresh=[false for _ in 1:length(MLFMMTrees.tree(receivetree).nodes)]
aggregationlist = _initialize_aggregationlist(sourcetree, transferlist)
disaggregationlist= _initialize_disaggregationlist(receivetree, transferlist)


verbose && @info "Allocate  node spectra"
# ensure that levelcutoffparameters are the same as for sourcestruct
levelcutoffparameters = _initializelevelcutoffparameters(receivetree, expectedaccuracy, wavenumber)
sourcelevels=MLFMMTrees.levels(sourcetree)
receivelevels=MLFMMTrees.levels(receivetree)
minlevels= minimum([length(sourcelevels), length(receivelevels)])
levelcutoffparameters[end - minlevels + 1 : end] = sourcestruct.levelcutoffparameters[end - minlevels + 1 : end] 
# nodespectra = _allocatenodepattern(PlaneWaveSpectrum{Complex{R}}, receivetree, levelcutoffparameters, transferlist)
nodespectra = _allocatenodepattern(PlaneWaveSpectrum{Complex{R}}, receivetree, levelcutoffparameters, minreceivelevel)
tmpmatrix=similar(nodespectra[end].Eθ)

verbose && @info "Assemble basis spectra"
leafnodeindices = MLFMMTrees.leafs(receivetree)
basisfunctionpatterns = _initializebasisfunctionweightingpatterns(receivetree, testfunctions, levelcutoffparameters, wavenumber, verbose)
bvector=Vector{Complex{R}}(undef, length(basisfunctionpatterns))


verbose && @info "Assemble interpolators"
levelinterpolators = _initializelevelinterpolators(R, receivetree, θstencilsize, ϕstencilsize, levelcutoffparameters; minlevel = minreceivelevel)
phaseshifttoparent = _initializephaseshifttoparent(R, receivetree,levelcutoffparameters, wavenumber; minlevel = minreceivelevel)

verbose && println("------------------------------")
return MLFMMReceive(receivetree, nodespectra, basisfunctionpatterns, levelcutoffparameters, levelinterpolators, phaseshifttoparent, transferlist, adjoint_transferlist, transferplan, nodeisfresh, leafnodeindices, minreceivelevel, minsourcetranslationlevel, receivetranslationnodes, firetranslationnodes, aggregationlist, disaggregationlist, bvector, verbose, tmpmatrix)
end

"""
    _initialize_tree(points::AbstractArray{SVector{3,R}}, minhalfsize::R) where{R<:Real}

Assemble 'MLFMMTree' based on point cloud and half length of leave boxes
"""
function _initialize_tree(points::AbstractArray{SVector{3,R}}, minhalfsize::R) where{R<:Real}
    rootcenter, rootsize = MLFMMTrees.getboundingbox(points)

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2=log2(rootsize/minhalfsize)
    halfsize_powerof2=maximum([0, Integer(ceil(halfsize_powerof2))])
    rootsize=2^halfsize_powerof2*minhalfsize

    root_center = SVector{3,Float64}(rootcenter)
    return ElectromagneticFieldRepresentations.MLFMMTrees.MLFMMTree(root_center, points, rootsize, minhalfsize)
end

"""
    _getpoints(basisfunctions)

Return list of `SVector{3}` representing the reference (i.e., center) locations of the `basisfunctions`
"""
function _getpoints(basisfunctions::AbstractSurfaceCurrentDensity)
    return BEAST.positions(functionspace(basisfunctions))
end
function _getpoints(basisfunctions::ElmagSurfaceCurrentDensity)
    pos=BEAST.positions(functionspace(basisfunctions))
    return [pos;pos]
end
function _getpoints(functionspace)
    return functionspace.points
end
function _getpoints(functionspace::BEAST.Space{R}) where{R<:Real}
    return BEAST.positions(functionspace)
end
function _getpoints(functionspace::AbstractArray{D}) where{D<:AbstractDipole}
    return [dipole.pos for dipole in functionspace][:]
end

"""
    _getsourcefunctions(basisfunctions)

Return list of individual source functions
"""
function _getsourcefunctions(basisfunctions::NamedTuple{(:sourcefunctions, :points), Tuple{V1,V2}}) where{V1,V2}
    return basisfunctions.sourcefunctions
end
function _getsourcefunctions(functionspace::AbstractArray{D}) where{D<:AbstractDipole}
    dipoles =deepcopy(functionspace)
    for (k, dipole) in enumerate functionspace
        dipoles[k] = D([0;0;0], dipole.dir, dipole.mag)
    end
end

"""
    _cutoffparameter(diameter, expectedaccuracy, k0)

Return cutoff parameter (i.e., maximum spherical mode index L) to represent sources inside sphere of given `diameter` with expected accuracy.
"""
function _cutoffparameter(diameter, expectedaccuracy, k0)
    kd=k0*diameter
    digitsofaccuracy=-log10(expectedaccuracy)
    L=ceil(Int, kd + 1.8 * digitsofaccuracy^(2 / 3) * (kd)^(1 / 3)) + 2
    return L
end

"""
    _initializelevelcutoffparameters(tree::MLFMMTrees.AbstractMLFMMTree, expectedaccuracy, wavenumber)

Return list of cutoff parameters for each level in `tree` to represent the node farfields with expected accuracy.
"""
function _initializelevelcutoffparameters(tree::MLFMMTrees.AbstractMLFMMTree, expectedaccuracy, wavenumber)
    levels=MLFMMTrees.levels(tree)
    cutoffparameters=Array{Int}(undef, length(levels))

    for level in levels
        node=MLFMMTrees.tree(tree).nodes[MLFMMTrees.nodesatlevel(tree,level)][1]
        halfsize=node.data.halfsize
        diameter=sqrt(3*halfsize^2)*2
        cutoffparameters[level]=_cutoffparameter(diameter, expectedaccuracy, wavenumber)
    end

    return cutoffparameters
end

# function _allocatenodefarfields(R::Type{<:Real}, tree::MLFMMTrees.AbstractMLFMMTree, cutoffparameters::Vector{<:Integer})
#     nodefarfields=Vector{FarfieldPattern{Complex{R}}}(undef, length(MLFMMTrees.tree(tree).nodes))
    
#     levels=MLFMMTrees.levels(tree)
#     for level in levels
#         L= cutoffparameters[level]
#         for node in MLFMMTrees.nodesatlevel(tree, level)
#             Eθ=Array{Complex{R}}(undef, L+1, 2L+2)
#             Eϕ=Array{Complex{R}}(undef, L+1, 2L+2)
#             nodefarfields[node]=FarfieldPattern(L, Eθ, Eϕ)
#         end
#     end

#     return nodefarfields
# end


"""
    _allocatenodepattern(P::Type{S}, tree::MLFMMTrees.AbstractMLFMMTree, cutoffparameters::Vector{<:Integer}, [transferlist::Vector{Vector{Int}}, minlevelel::Int]) where {R<:Real, S<:Union{PlaneWaveSpectrum{Complex{R}}, FarfieldPattern{Complex{R}}}}

Allocate memory for intermediate representation of far field patterns for each node in tree. 
    
If `transferlist` is provided, only node patterns which are relevant for translations are stored.
If `minlevel` is provided, only node patterns up to minlevel are stored.
"""
function _allocatenodepattern(P::Type{S}, tree::MLFMMTrees.AbstractMLFMMTree, cutoffparameters::Vector{<:Integer}, transferlist::Vector{Vector{Int}}) where {R<:Real, S<:Union{PlaneWaveSpectrum{Complex{R}}, FarfieldPattern{Complex{R}}}}
# Only allocate patterns for receiving nodes and their children
    nodepattern=Vector{P}(undef, length(MLFMMTrees.tree(tree).nodes))
    
    levels=MLFMMTrees.levels(tree)
    for level in levels
        L= cutoffparameters[level]
        for node in MLFMMTrees.nodesatlevel(tree, level)
            if isempty(transferlist[node])
                nodepattern[node]=P(0, Array{Complex{R}}(undef,0,0), Array{Complex{R}}(undef,0,0))
            else
                Eθ=Array{Complex{R}}(undef, L+1, 2L+2)
                Eϕ=Array{Complex{R}}(undef, L+1, 2L+2)
                nodepattern[node]=P(L, Eθ, Eϕ)
            end
        end
    end
    for node::Int in eachindex(nodepattern)
        if nodepattern[node].L != 0
            for childnode::Int in MLFMMTrees.DepthFirstIterator(tree, node)
                if (nodepattern[childnode].L == 0)
                    L= cutoffparameters[MLFMMTrees.level(tree, childnode)]
                    Eθ=Array{Complex{R}}(undef, L+1, 2L+2)
                    Eϕ=Array{Complex{R}}(undef, L+1, 2L+2)
                    nodepattern[childnode]=P(L, Eθ, Eϕ)
                end
            end
        end
    end

    return nodepattern
end
function _allocatenodepattern(P::Type{S}, tree::MLFMMTrees.AbstractMLFMMTree, cutoffparameters::Vector{<:Integer}) where {R<:Real, S<:Union{PlaneWaveSpectrum{Complex{R}}, FarfieldPattern{Complex{R}}}}
    nodepattern=Vector{P}(undef, length(MLFMMTrees.tree(tree).nodes))
    
    levels=MLFMMTrees.levels(tree)
    for level in levels
        L= cutoffparameters[level]
        for node in MLFMMTrees.nodesatlevel(tree, level)
            Eθ=Array{Complex{R}}(undef, L+1, 2L+2)
            Eϕ=Array{Complex{R}}(undef, L+1, 2L+2)
            nodepattern[node]=P(L, Eθ, Eϕ)
        end
    end

    return nodepattern
end
function _allocatenodepattern(P::Type{S}, tree::MLFMMTrees.AbstractMLFMMTree, cutoffparameters::Vector{<:Integer}, minlevel::Integer) where {R<:Real, S<:Union{PlaneWaveSpectrum{Complex{R}}, FarfieldPattern{Complex{R}}}}
    nodepattern=Vector{P}(undef, length(MLFMMTrees.tree(tree).nodes))
    
    levels=MLFMMTrees.levels(tree)
    for level in minlevel: length(levels)
        L= cutoffparameters[level]
        for node in MLFMMTrees.nodesatlevel(tree, level)
            Eθ=Array{Complex{R}}(undef, L+1, 2L+2)
            Eϕ=Array{Complex{R}}(undef, L+1, 2L+2)
            nodepattern[node]=P(L, Eθ, Eϕ)
        end
    end

    return nodepattern
end

"""
    _initializelevelinterpolators(R::Type{<:Real}, tree::MLFMMTrees.AbstractMLFMMTree, orderθ::Integer, orderϕ::Integer,cutoffparameters::Vector{<:Integer}; minlevel::Int=0)

Return list of interpolators to aggregate from each level to the respective parent level.
"""
function _initializelevelinterpolators(R::Type{<:Real}, tree::MLFMMTrees.AbstractMLFMMTree, orderθ::Integer, orderϕ::Integer,cutoffparameters::Vector{<:Integer}; minlevel::Int=0)
    levels=MLFMMTrees.levels(tree)

    interplatortype=ParallelPlannedLocalInterpolation{orderθ, orderϕ, R}
    interplatortype2=ParallelPlannedLocalInterpolation{R}

    levelinterpolators=Vector{interplatortype}(undef, length(levels)-1)
    for level in levels[1:end-1]
        level< minlevel && continue

        Lold=cutoffparameters[level+1]
        Lnew=cutoffparameters[level]
        levelinterpolators[level]=interplatortype2(Lold,Lnew, orderθ=orderθ, orderϕ=orderϕ)
        # levelinterpolators[level]=PlannedLocalThetaGlobalPhiInterpolation{R}(Lold,Lnew, orderθ=minimum([orderθ,Lold]))
        # levelinterpolators[level]=GlobalInterpolation{R}(Lold,Lnew) 
        # if Lold >= maximum([orderθ; orderϕ])
        #     levelinterpolators[level]=PlannedLocalInterpolation{R}(Lold,Lnew, orderθ=orderθ, orderϕ=orderϕ)
        # else
        #     levelinterpolators[level]=GlobalInterpolation{R}(Lold,Lnew) 
        # end      
    end
    return levelinterpolators
end
#TODO: only initialize required phaseshifts and share between source and receive
function _initializephaseshifttoparent(T::Type{<:Real}, tree::MLFMMTrees.AbstractMLFMMTree,cutoffparameters::Vector{<:Integer}, k0::Real; minlevel::Int=0)
    levels=MLFMMTrees.levels(tree)
    phaseshifttoparent=Array{Matrix{Complex{T}}}(undef, 8 , length(levels))
    for level in levels[2:end]
        level < minlevel && continue

        L=cutoffparameters[level-1]
        # phaseshiftmatrix=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        node=MLFMMTrees.tree(tree).nodes[MLFMMTrees.nodesatlevel(tree,level)[1]]
        halfsize=node.data.halfsize

        R=[halfsize, halfsize, halfsize]
        phaseshifttoparent[1,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[1,level], R, k0)        

        R=[-halfsize, halfsize, halfsize]
        phaseshifttoparent[2,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[2,level], R, k0)

        R=[halfsize,-halfsize, halfsize]
        phaseshifttoparent[3,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[3,level], R, k0)

        R=[-halfsize,-halfsize, halfsize]
        phaseshifttoparent[4,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[4,level], R, k0)

        R=[halfsize,halfsize, -halfsize]
        phaseshifttoparent[5,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[5,level], R, k0)

        R=[-halfsize, halfsize, -halfsize]
        phaseshifttoparent[6,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[6,level], R, k0)

        R=[halfsize,-halfsize, -halfsize]
        phaseshifttoparent[7,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[7,level], R, k0)

        R=[-halfsize,-halfsize, -halfsize]
        phaseshifttoparent[8,level]=Matrix{Complex{T}}(undef, L + 1, 2 * L + 2)
        _phaseshiftmatrix!(phaseshifttoparent[8,level], R, k0)
    end
    return phaseshifttoparent
end

struct typeindicator{T}
end
#TODO: share memory beween levels
"""
    _initializebasisfunctionfarfields(tree::MLFMMTrees.AbstractMLFMMTree, basisfunctions::S, cutoffparameters::Vector{I}, k0::R, verbose::Bool) where {R<:Real, S<:AbstractSurfaceCurrentDensity, I<:Int}

Return list of far fields for every basis function on leaf level.
"""
function _initializebasisfunctionfarfields(tree::MLFMMTrees.AbstractMLFMMTree, basisfunctions::S, cutoffparameters::Vector{I}, k0::R, verbose::Bool) where {R<:Real, S<:AbstractSurfaceCurrentDensity, I<:Int}
    basisfunctionfarfields=Vector{FarfieldPattern{Complex{R}}}(undef, numfunctions(basisfunctions))
    L= cutoffparameters[end]
    _, θvec, ϕvec = samplingrule(L)
    nθ = length(θvec)
    nϕ = length(ϕvec)

    pts = [point(cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)) for ϕ in ϕvec for θ in θvec]

    farfieldmatrices=individualfarfields(basisfunctions, pts, k0, quadstrat=BEAST.SingleNumQStrat(1), verbose=false)::Matrix{SVector{3,Complex{R}}}
    cost = cos.(θvec)
    sint = sin.(θvec)

    cosp = cos.(ϕvec)
    sinp = sin.(ϕvec)
    
    ff = Matrix{SVector{3,Complex{R}}}(undef, nθ, nϕ)
    for functionindex in eachindex(basisfunctionfarfields)
        ff .= reshape(view(farfieldmatrices,:,functionindex), nθ, nϕ)

        basisfunctionfarfields[functionindex]=FarfieldPattern{Complex{R}}(L, Matrix{Complex{R}}(undef, length(θvec), length(ϕvec)), Matrix{Complex{R}}(undef, length(θvec), length(ϕvec)))
        Eθ = basisfunctionfarfields[functionindex].Eθ
        Eϕ = basisfunctionfarfields[functionindex].Eϕ
        for kθ in eachindex(θvec), kϕ in eachindex(ϕvec)
        eθ = SVector{3,R}(cost[kθ] * cosp[kϕ], cost[kθ] * sinp[kϕ], -sint[kθ])
        eϕ = SVector{3,R}(-sinp[kϕ], cosp[kϕ], zero(eltype(θvec)))

        Eθ[kθ, kϕ] = Complex{R}(udot(eθ, ff[kθ, kϕ]))
        Eϕ[kθ, kϕ] = Complex{R}(udot(eϕ, ff[kθ, kϕ]))  
        end
            
    end
    phaseshiftmatrix=Array{Complex{R}}(undef, nθ, nϕ)
        for leafnode::Int in MLFMMTrees.leafs(tree)
            _phaseshiftmatrix!(phaseshiftmatrix, MLFMMTrees.center(tree, leafnode), k0)
            for functionindex::Int in tree(leafnode).data.values
                basisfunctionfarfields[functionindex].Eθ .*= (phaseshiftmatrix)
                basisfunctionfarfields[functionindex].Eϕ .*= (phaseshiftmatrix)
            end
        end

    return basisfunctionfarfields
end
function _initializebasisfunctionfarfields(tree::TT, basisfunctions,
    cutoffparameters::Vector{I}, k0::R, verbose::Bool) where{R<:Real, I<:Integer, TT<:MLFMMTrees.AbstractMLFMMTree}

    sourcefunctions=_getsourcefunctions(basisfunctions)
    points=_getpoints(basisfunctions)

    basisfunctionfarfields=Vector{FarfieldPattern{Complex{R}}}(undef, length(sourcefunctions)) 
    L= cutoffparameters[end]
    _, θvec, ϕvec = samplingrule(L)
    nθ = length(θvec)
    nϕ = length(ϕvec)
    phaseshiftmatrix=Array{Complex{R}}(undef, nθ, nϕ)
    
    for (f_index, sourcefunction) in enumerate(sourcefunctions)
        basisfunctionfarfields[f_index]= _convertrepresentation_and_resample(FarfieldPattern{Complex{R}}, sourcefunction, L, k0)
        _phaseshiftmatrix!(phaseshiftmatrix, -points[f_index], k0)
        basisfunctionfarfields[f_index].Eθ .*= (phaseshiftmatrix)
        basisfunctionfarfields[f_index].Eϕ .*= (phaseshiftmatrix)

    end
   
    for leafnode::Int in MLFMMTrees.leafs(tree)
        _phaseshiftmatrix!(phaseshiftmatrix, MLFMMTrees.center(tree, leafnode), k0)
        for functionindex::Int in tree(leafnode).data.values
            basisfunctionfarfields[functionindex].Eθ .*= (phaseshiftmatrix)
            basisfunctionfarfields[functionindex].Eϕ .*= (phaseshiftmatrix)
        end
    end

    return basisfunctionfarfields
end

"""
    _initializebasisfunctionweightingpatterns(tree::MLFMMTrees.AbstractMLFMMTree, functionspace,
    cutoffparameters::Vector{<:Integer}, k0::R, verbose::Bool) where{R<:Real}

Return list of far fields for every test function on leaf level stored in reverted direction feasible for testing incident plane wave spectra. 
"""
function _initializebasisfunctionweightingpatterns(tree::MLFMMTrees.AbstractMLFMMTree, functionspace,
    cutoffparameters::Vector{<:Integer}, k0::R, verbose::Bool) where{R<:Real}
    basisfunctionfarfields=_initializebasisfunctionfarfields(tree, functionspace, cutoffparameters, -k0, verbose)
    for (index, farfield) in enumerate(basisfunctionfarfields)
        # basisfunctionfarfields[index]=revertdirection(farfield)
        basisfunctionfarfields[index]=(farfield)
    end
    return basisfunctionfarfields
end

"""
    _aggregate_leafnodes!(A::MLFMMSource)

Fill storage for leaf node patterns due to excitation vector `A.xvector`.
"""
function _aggregate_leafnodes!(A::MLFMMSource)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)

    # for leafnode in A.leafnodeindices
    Threads.@threads for leafnode in A.leafnodeindices
        reset =true
        for functionindex::Int in tree(leafnode).data.values::Vector{Int}
            _muladd!(A.nodefarfields[leafnode],A.basisfunctionfarfields[functionindex], A.xvector[functionindex], reset=reset)
            reset=false
        end         
    end
end
#TODO: store children per node to remove dependency on ClusteTrees.children 
"""
    _aggregate_children!(A::MLFMMSource, parentnode)

Fill pattern storage of parentnode with aggregated pattern from all its children.
"""
function _aggregate_children!(A::MLFMMSource, parentnode)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)

    level = MLFMMTrees.level(tree, parentnode)

            interpolator= A.levelinterpolators[level]
            reset =true
            
            for child::Int in ClusterTrees.children(tree, parentnode)
                
                sector= tree.nodes[child].data.sector+1

                _interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[child].Eθ, interpolator)
                _muladd!(A.nodefarfields[parentnode].Eθ, interpolator.finalstorage[Threads.threadid()], (A.phaseshifttoparent[sector,level+1]), reset=reset)

                _interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[child].Eϕ, interpolator)
                _muladd!(A.nodefarfields[parentnode].Eϕ, interpolator.finalstorage[Threads.threadid()], (A.phaseshifttoparent[sector,level+1]), reset=reset)

                reset = false
            end
end

"""
_aggregate_to_minlevel!(A::MLFMMSource, [x::AbstractVector]; min_aggregationlevel::Integer=0)

Aggregate `A` up to min_aggregationlevel. 
"""
function _aggregate_to_minlevel!(A::MLFMMSource, x; min_aggregationlevel::Integer=0)
    A.xvector .= x
    _aggregate_to_minlevel!(A, min_aggregationlevel = min_aggregationlevel)
end
function _aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer=0)
    A.verbose && @info "Aggregate node far fields"
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)

    _aggregate_leafnodes!(A)
    
    levels=MLFMMTrees.levels(tree)
    for level in reverse(maximum([min_aggregationlevel,1]):length(levels))
        for parentnode::Int in MLFMMTrees.nodesatlevel(tree, level)
            MLFMMTrees.isleaf(tree, parentnode) && continue            
            _aggregate_children!(A, parentnode)
        end
    end

end

"""
    _aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})

Aggregate `A` according to aggregationlist.

The aggregationlist specifies for each level which nodes shall be aggregated.
"""
function _aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})
    A.verbose && @info "Aggregate node far fields"

    _aggregate_leafnodes!(A)
    
    for level in reverse(eachindex(aggregationlist))
        # for node in aggregationlist[level]
        Threads.@threads for node in aggregationlist[level]
            _aggregate_children!(A, node)
        end
    end

end

"""
    _adjoint_aggregate_leafnodes!(A::MLFMMSource)

Perform the adjoint operation (i.e., conjugate of transposed operation) to `_aggregate_leafnodes!`
"""
function _adjoint_aggregate_leafnodes!(A::MLFMMSource)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)


    # for leafnode in A.leafnodeindices
    #     pws=A.nodefarfields[leafnode]
    #     for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
    #         ff=A.basisfunctionfarfields[basisfunctionindex]
    #         # receivestruct.bvector[probeindex] = _ewaldintegral(ff, pws, w, L)
    #         A.tmpmatrix .= (pws.Eθ .* conj.(ff.Eθ)) .+ (pws.Eϕ .* conj.(ff.Eϕ)) 
    #         A.xvector[basisfunctionindex] = sum(A.tmpmatrix) 
    #     end
    # end
    Threads.@threads for leafnode in A.leafnodeindices
        pws=A.nodefarfields[leafnode]
        # pwsEθ = reshape(pws.Eθ,:)
        # pwsEϕ = reshape(pws.Eϕ,:)
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff=A.basisfunctionfarfields[basisfunctionindex]
            # ffEθ = reshape(ff.Eθ,:)
            # ffEϕ = reshape(ff.Eϕ,:)
            A.xvector[basisfunctionindex] = dot(ff.Eθ, pws.Eθ) .+  dot(ff.Eϕ, pws.Eϕ)
        end
    end
end


"""
    _transpose_aggregate_leafnodes!(A::MLFMMSource)

Perform the transpose operation to `_aggregate_leafnodes!`
"""
function _transpose_aggregate_leafnodes!(A::MLFMMSource)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)


    # for leafnode in A.leafnodeindices
    #     pws=A.nodefarfields[leafnode]
    #     for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
    #         ff=A.basisfunctionfarfields[basisfunctionindex]
    #         # receivestruct.bvector[probeindex] = _ewaldintegral(ff, pws, w, L)
    #         A.tmpmatrix .= (pws.Eθ .* conj.(ff.Eθ)) .+ (pws.Eϕ .* conj.(ff.Eϕ)) 
    #         A.xvector[basisfunctionindex] = sum(A.tmpmatrix) 
    #     end
    # end
    Threads.@threads for leafnode in A.leafnodeindices
        pws=A.nodefarfields[leafnode]
        # pwsEθ = reshape(pws.Eθ,:)
        # pwsEϕ = reshape(pws.Eϕ,:)
        for basisfunctionindex::Int in tree(leafnode).data.values::Vector{Int}
            ff=A.basisfunctionfarfields[basisfunctionindex]
            # ffEθ = reshape(ff.Eθ,:)
            # ffEϕ = reshape(ff.Eϕ,:)
            A.xvector[basisfunctionindex] = udot(ff.Eθ, pws.Eθ) .+  udot(ff.Eϕ, pws.Eϕ)
        end
    end
end

#TODO: store children per node to remove dependency on ClusteTrees.children 
"""
    _transpose_aggregate_children!(A, parentnode)

Perform transposed operation to `_aggregate_children!`
"""
function _transpose_aggregate_children!(A, parentnode)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}
    tree::treetype=MLFMMTrees.tree(A.tree)
    
    level = MLFMMTrees.level(tree, parentnode)
    interpolator= A.levelinterpolators[level]

            for child::Int in ClusterTrees.children(tree, parentnode)

                sector= tree.nodes[child].data.sector+1

                _muladd!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[parentnode].Eθ, A.phaseshifttoparent[sector,level+1], reset=true)
                _transpose_interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[child].Eθ, interpolator, reset= !A.nodeisfresh[child])

                _muladd!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[parentnode].Eϕ, A.phaseshifttoparent[sector,level+1], reset=true)
                _transpose_interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodefarfields[child].Eϕ, interpolator, reset= !A.nodeisfresh[child])
            
                A.nodeisfresh[child]=true
            end
end

#TODO: avoid complex conjugation for phaseshifts (use different sectors). 
"""
    _adjoint_aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})

Perform adjoint operation (i.e., complex conjugate of transposed operation) of ` _adjoint_aggregate!`
"""
function _adjoint_aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})
    A.verbose && @info "Adjoint aggregate node far fields"
    
    for level in eachindex(aggregationlist)
        aggregationlist[level]==[] && continue
        # for sector in 1:8
        Threads.@threads for sector in 1:8
            conj!(A.phaseshifttoparent[sector,level+1])
        end

        # for node in aggregationlist[level]
        Threads.@threads for node in aggregationlist[level]
            _transpose_aggregate_children!(A, node)
        end

        # for sector in 1:8
        Threads.@threads for sector in 1:8
            conj!(A.phaseshifttoparent[sector,level+1])
        end
    end

    _adjoint_aggregate_leafnodes!(A)

end

"""
    _transpose_aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})

Perform transpose operation of `_aggregate!`
"""
function _transpose_aggregate!(A::MLFMMSource, aggregationlist::Vector{Vector{Int}})
    A.verbose && @info "Adjoint aggregate node far fields"
    
    for level in eachindex(aggregationlist)
        aggregationlist[level]==[] && continue

        # for node in aggregationlist[level]
        Threads.@threads for node in aggregationlist[level]
            _transpose_aggregate_children!(A, node)
        end
    end

    _transpose_aggregate_leafnodes!(A)

end

"""
    _disaggregate_leafnodes!(receivestruct::MLFMMReceive)

Test leaf node patterns with all test functions of corresponding leafnode and store result in `receivestruct.bvector`
"""
function _disaggregate_leafnodes!(receivestruct::MLFMMReceive)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}    
    receivetree::treetype = MLFMMTrees.tree(receivestruct.tree)

    Threads.@threads for leafnode in receivestruct.leafnodeindices
        pws=receivestruct.nodespectra[leafnode]

        for probeindex::Int in receivetree(leafnode).data.values::Vector{Int}
            ff=receivestruct.basisfunctionpatterns[probeindex]

            receivestruct.bvector[probeindex] = udot(pws.Eθ, ff.Eθ) .+  udot(pws.Eϕ, ff.Eϕ)
        end
    end

end

"""
    _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)

Perform ajoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate_leafnodes!`
"""
function _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}    
    receivetree::treetype = MLFMMTrees.tree(A.tree)

    Threads.@threads for leafnode in A.leafnodeindices
        reset=true
        for functionindex::Int in receivetree(leafnode).data.values::Vector{Int}
            if !reset
                A.nodespectra[leafnode].Eθ .+= A.bvector[functionindex] .* conj.(A.basisfunctionpatterns[functionindex].Eθ)
                A.nodespectra[leafnode].Eϕ .+= A.bvector[functionindex] .* conj.(A.basisfunctionpatterns[functionindex].Eϕ) 
            else
                A.nodespectra[leafnode].Eθ .= A.bvector[functionindex] .* conj.(A.basisfunctionpatterns[functionindex].Eθ)
                A.nodespectra[leafnode].Eϕ .= A.bvector[functionindex] .* conj.(A.basisfunctionpatterns[functionindex].Eϕ) 
                reset=false 
            end
        end         
    end
    
end

"""
    _transpose_disaggregate_leafnodes!(A::MLFMMReceive)

Perform transposed operation to `disaggregate_leafnodes!`
"""
function _transpose_disaggregate_leafnodes!(A::MLFMMReceive)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}    
    receivetree::treetype = MLFMMTrees.tree(A.tree)

    Threads.@threads for leafnode in A.leafnodeindices
        reset=true
        for functionindex::Int in receivetree(leafnode).data.values::Vector{Int}
            _muladd!(A.nodespectra[leafnode], A.basisfunctionpatterns[functionindex], A.bvector[functionindex], reset=reset)
            reset=false
        end         
    end
    
end

#TODO: store children per node in a list to remove dependency on ClusteTrees.children 
"""
    _disaggregate_children!(receivestruct, parentnode)

Disaggregate pattern of parentnode to all its children and store the resulting patterns in `receivestruct.nodespectra[child]`
"""
function _disaggregate_children!(receivestruct, parentnode)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}    
    receivetree::treetype = MLFMMTrees.tree(receivestruct.tree)

    level = MLFMMTrees.level(receivetree, parentnode)
    interpolator= receivestruct.levelinterpolators[level]
    for child::Int in ClusterTrees.children(receivetree, parentnode)
        sector= receivetree.nodes[child].data.sector+1

        _muladd!(interpolator.finalstorage[Threads.threadid()], receivestruct.nodespectra[parentnode].Eθ, receivestruct.phaseshifttoparent[sector,level+1], reset=true)
        _adjoint_interpolatematrix!(interpolator.finalstorage[Threads.threadid()], receivestruct.nodespectra[child].Eθ, interpolator, reset= !receivestruct.nodeisfresh[child])

        _muladd!(interpolator.finalstorage[Threads.threadid()], receivestruct.nodespectra[parentnode].Eϕ, receivestruct.phaseshifttoparent[sector,level+1], reset=true)
        _adjoint_interpolatematrix!(interpolator.finalstorage[Threads.threadid()], receivestruct.nodespectra[child].Eϕ, interpolator, reset= !receivestruct.nodeisfresh[child])
        
        receivestruct.nodeisfresh[child] =true
    end
end

#TODO: store children per node to remove dependency on ClusteTrees.children   
"""
    _transpose_disaggregate_children!(A, parentnode)

Perform transposed operation to `disaggrgate_children!`
"""
function _transpose_disaggregate_children!(A, parentnode)
    treetype = MLFMMTrees.PBMLFMMTree{ElectromagneticFieldRepresentations.MLFMMTrees.Data{3, Float64}, 3, Float64}    
    receivetree::treetype = MLFMMTrees.tree(A.tree)
    level = MLFMMTrees.level(receivetree, parentnode)

    interpolator= A.levelinterpolators[level]
            reset =true
                for child::Int in ClusterTrees.children(receivetree, parentnode)
                    
                    sector= receivetree.nodes[child].data.sector+1

                    _interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodespectra[child].Eθ, interpolator)
                    _muladd!(A.nodespectra[parentnode].Eθ, interpolator.finalstorage[Threads.threadid()], (A.phaseshifttoparent[sector,level+1]), reset=reset)

                    _interpolatematrix!(interpolator.finalstorage[Threads.threadid()], A.nodespectra[child].Eϕ, interpolator)
                    _muladd!(A.nodespectra[parentnode].Eϕ, interpolator.finalstorage[Threads.threadid()], (A.phaseshifttoparent[sector,level+1]), reset=reset)

                    reset = false
                end
end

"""
    _disaggregate!(receivestruct::MLFMMReceive)

Perform all disaggregations according to `receivestruct.disaggregationlist`.
The `receivestruct.disaggregationlist` stores all nodes at each level which shall perform a disaggregation.
"""
function _disaggregate!(receivestruct::MLFMMReceive)
    receivestruct.verbose && @info "Disaggregate node spectra "

    
    for level in eachindex(receivestruct.disaggregationlist)
        receivestruct.disaggregationlist[level]==[] && continue
        # for sector in 1:8
        Threads.@threads for sector in 1:8
            conj!(receivestruct.phaseshifttoparent[sector,level+1])
        end

        # for node in receivestruct.disaggregationlist[level]
        Threads.@threads for node in receivestruct.disaggregationlist[level]
            _disaggregate_children!(receivestruct, node)
        end

        # for sector in 1:8
        Threads.@threads for sector in 1:8
            conj!(receivestruct.phaseshifttoparent[sector,level+1])
        end

    end

    _disaggregate_leafnodes!(receivestruct)

end

"""
    _adjoint_disaggregate!(A::MLFMMReceive, [y::AbstractVector])

Perform adjoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate!`#

If no vector `y` is given, the content in `A.bvector` is used for adjoint disaggregation. 
Otherwise, `A.bvector` is overwritten by `y` before adjoint disaggregation.
"""
function _adjoint_disaggregate!(A::MLFMMReceive, y::AbstractVector)
    A.bvector .= y
    _adjoint_disaggregate!(A)
end
function _adjoint_disaggregate!(A::MLFMMReceive)
    A.verbose && @info "Adjoint disaggregate node spectra "
    
    _adjoint_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))

        # for node in A.disaggregationlist[level]
        Threads.@threads for node in A.disaggregationlist[level]
            # _adjoint_disaggregate_children!(A, node)
            _transpose_disaggregate_children!(A, node)
        end

    end
    # _adjoint_disaggregate_leafnodes!(A)

end

"""
    _transpose_disaggregate!(A::MLFMMReceive, [y::AbstractVector])

    Perform transposed operation of `_disaggregate!`

If no vector `y` is given, the content in `A.bvector` is used for transposed disaggregation. 
Otherwise, `A.bvector` is overwritten by `y` before transposed disaggregation.
"""
function _transpose_disaggregate!(A::MLFMMReceive, y::AbstractVector)
    A.bvector .= y
    _transpose_disaggregate!(A)
end
function _transpose_disaggregate!(A::MLFMMReceive)
    A.verbose && @info "Transpose disaggregate node spectra "
    
    _transpose_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))
        A.disaggregationlist[level]==[] && continue
        
        Threads.@threads for sector in 1:8
            conj!(A.phaseshifttoparent[sector,level+1])
        end

        
        Threads.@threads for node in A.disaggregationlist[level]
            _transpose_disaggregate_children!(A, node)
        end

        Threads.@threads for sector in 1:8
            conj!(A.phaseshifttoparent[sector,level+1])
        end

    end

end

"""
    _initialize_aggregationlist(TXtree::MLFMMTrees.AbstractMLFMMTree, transferlist)

Return a list of nodes per level which shall be aggregated to be able to perform all transfers provided in `transferlist`.
"""
function _initialize_aggregationlist(TXtree::MLFMMTrees.AbstractMLFMMTree, transferlist)
    sourcetree = MLFMMTrees.tree(TXtree)
    numlevels=length(sourcetree.nodesatlevel)
    # aggregationlist = Vector{Vector{Int}}(undef, numlevels)
    aggregationlist = [Vector{Int}([]) for k in 1:numlevels]
    for transmitlist in transferlist
        for translatenode in transmitlist
            for childnode in MLFMMTrees.DepthFirstIterator(sourcetree, translatenode)
                MLFMMTrees.isleaf(sourcetree, childnode) && continue
                level = MLFMMTrees.level(sourcetree, childnode)
                push!(aggregationlist[level], childnode)
            end
        end
    end

    for level in eachindex(aggregationlist)
        aggregationlist[level]=sort!(unique(aggregationlist[level]))
    end
    return aggregationlist
end


"""
    _initialize_disaggregationlist(RXtree::MLFMMTrees.AbstractMLFMMTree, transferlist)

Return a list of nodes per level which shall be disaggregated to be able to perform all  disaggregations from nodes which receive transfers provided in `transferlist`.
"""
function _initialize_disaggregationlist(RXtree::MLFMMTrees.AbstractMLFMMTree, transferlist)
    receivetree = MLFMMTrees.tree(RXtree)
    numlevels=length(receivetree.nodesatlevel)
    # disaggregationlist = Vector{Vector{Int}}(undef, numlevels)
    disaggregationlist = [Vector{Int}([]) for k in 1:numlevels]
    for receivenode in eachindex(transferlist)
        if transferlist[receivenode] != []
            for childnode in MLFMMTrees.DepthFirstIterator(receivetree, receivenode)
                MLFMMTrees.isleaf(receivetree, childnode) && continue
                level = MLFMMTrees.level(receivetree, childnode)
                push!(disaggregationlist[level], childnode)
            end
        end
    end
    for level in eachindex(disaggregationlist)
        disaggregationlist[level]=sort!(unique(disaggregationlist[level]))
    end
    return disaggregationlist

end

"""
    _initialize_transferlist(TXtree::MLFMMTrees.AbstractMLFMMTree, RXtree::MLFMMTrees.AbstractMLFMMTree; num_bufferboxes::Integer=1)

For each node in `RXtree`, return a list of all transmit nodes in `TXtree` which shall perform a transfer to the receive node.
"""
function _initialize_transferlist(TXtree::MLFMMTrees.AbstractMLFMMTree, RXtree::MLFMMTrees.AbstractMLFMMTree; num_bufferboxes::Integer=1)
    receivetree = MLFMMTrees.tree(RXtree)
    sourcetree = MLFMMTrees.tree(TXtree)

    num_sourcelevels= length(sourcetree.nodesatlevel)
    num_receivelevels= length(receivetree.nodesatlevel)
    diff_levels=num_receivelevels-num_sourcelevels

    transferlist=[Int[] for _ in 1:MLFMMTrees.numberofnodes(receivetree)]
    adjoint_transferlist=[Int[] for _ in 1:MLFMMTrees.numberofnodes(sourcetree)]
    minreceivelevel::Int= typemax(Int)
    minsourcetranslationlevel::Int= typemax(Int)

    for receivenode in MLFMMTrees.DepthFirstIterator(receivetree, MLFMMTrees.root(receivetree))
        equivalent_sourcelevel = MLFMMTrees.level(receivetree, receivenode) - diff_levels
        equivalent_sourcelevel < 1 && continue


        for sourcenode in sourcetree.nodesatlevel[equivalent_sourcelevel]

            if MLFMMTrees.transferMLFMM(sourcetree,  receivetree, sourcenode, receivenode, num_bufferboxes)
                append!(transferlist[receivenode], sourcenode)
                append!(adjoint_transferlist[sourcenode], receivenode)
                minsourcetranslationlevel=minimum([minsourcetranslationlevel, MLFMMTrees.level(sourcetree, sourcenode)-1])
                minreceivelevel=minimum([minreceivelevel, MLFMMTrees.level(receivetree, receivenode)-1])
            end
        end
    end
    receivetranslationnodes =[k for k in eachindex(transferlist) if transferlist[k] != []]
    firetranslationnodes =[k for k in eachindex(adjoint_transferlist) if adjoint_transferlist[k] != []]
    return transferlist, adjoint_transferlist, minreceivelevel, minsourcetranslationlevel, receivetranslationnodes, firetranslationnodes
end

"""
    _initializetransferplan(transferlist::Vector{Vector{Int}}, sourcestruct::MLFMMSource, RXtree::MLFMMTrees.AbstractMLFMMTree;  ) where{P<:Type{<:AbstractTransfer}}

For each planned transfer, initialize and store a transferobject which holds all rewuired information to calculate the corresponding far field transfer.

## Keywordarguments:
transfertype::P=PlannedTransfer{Complex{R}} : specify type of the initialized tranfer object
"""
function _initializetransferplan(transferlist::Vector{Vector{Int}}, sourcestruct::MLFMMSource, RXtree::MLFMMTrees.AbstractMLFMMTree; transfertype::P=PlannedTransfer{Complex{R}} ) where{P<:Type{<:AbstractTransfer}}
    receivetree = MLFMMTrees.tree(RXtree)
    sourcetree = MLFMMTrees.tree(sourcestruct.tree)

    transferplan = [Vector{transfertype}(undef, length(sourcetree.nodes)) for _ in eachindex(transferlist)]

    Pℓstorage=[Vector{typeof(sourcestruct.wavenumber)}(undef, sourcestruct.levelcutoffparameters[1]+1) for _ in 1: Threads.nthreads()]

    Threads.@threads for receivenode in eachindex(transferlist)
        transfers= transferlist[receivenode]

        for transfernode in transfers
            Rvec = MLFMMTrees.center(receivetree, receivenode) - MLFMMTrees.center(sourcetree, transfernode)
            transferplan[receivenode][transfernode] =_initialize_transfer!(transfertype, Pℓstorage[Threads.threadid()], Rvec, sourcestruct.wavenumber, sourcestruct.levelcutoffparameters[MLFMMTrees.level(sourcetree, transfernode)])
         end

    end
    return transferplan
end

"""
    _transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

Perform all required transfers between the `sourcestruct` and the `receivestruct`.
"""
function _transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)
    transferlist = receivestruct.transferlist
    transferplan= receivestruct.transferplan
    sourcestruct.verbose && @info "Transfer source → receive"
    receivestruct.nodeisfresh .= false

    Threads.@threads for receivenode in receivestruct.receivetranslationnodes
        transfers=transferlist[receivenode]
            for transfernode in transfers
                translate!(receivestruct.nodespectra[receivenode], sourcestruct.nodefarfields[transfernode], transferplan[receivenode][transfernode], reset=!receivestruct.nodeisfresh[receivenode] )
                receivestruct.nodeisfresh[receivenode]=true
            end
    end
    nothing
end

"""
    _adjoint_transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

Perform the adjoint operator (i.e., complex conjugate of transposed operator) of `transfer!`
"""
function _adjoint_transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

    adjoint_transferlist = receivestruct.adjoint_transferlist
    transferplan= receivestruct.transferplan
    sourcestruct.verbose && @info "Adjoint transfer source ← receive"

    sourcestruct.nodeisfresh .= false

    Threads.@threads for transfernode in receivestruct.firetranslationnodes
        transfers=adjoint_transferlist[transfernode]
            for receivenode in transfers
                _adjoint_translate!(receivestruct.nodespectra[receivenode], sourcestruct.nodefarfields[transfernode], transferplan[receivenode][transfernode], reset=!sourcestruct.nodeisfresh[transfernode])
                sourcestruct.nodeisfresh[transfernode]=true
            end
    end
    nothing

end

"""
    _transpose_transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

Perform the transposed operator to `transfer!`
"""
function _transpose_transfer!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

    adjoint_transferlist = receivestruct.adjoint_transferlist
    transferplan= receivestruct.transferplan
    sourcestruct.verbose && @info "transpose transfer source ← receive"

    sourcestruct.nodeisfresh .= false

    Threads.@threads for transfernode in receivestruct.firetranslationnodes
        transfers=adjoint_transferlist[transfernode]
            for receivenode in transfers
                _transpose_translate!(receivestruct.nodespectra[receivenode], sourcestruct.nodefarfields[transfernode], transferplan[receivenode][transfernode], reset=!sourcestruct.nodeisfresh[transfernode])
                sourcestruct.nodeisfresh[transfernode]=true
            end
    end
    nothing

end

"""
    _forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive,[x])

Store the result of the matrix vector product in `receivestruct.bvector`

If no vector `x` is given, `sourcestruct.xvector` is used as excitation vector.
Otherwise, `sourcestruct.xvector` is overwritten by the content of `x` before the operation is performed
"""
function _forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)
    (sourcestruct.verbose || receivestruct.verbose) && @info "-------------------------\n   Evaluate forward operator\n-------------------------------"

    # _aggregate_to_minlevel!(sourcestruct; min_aggregationlevel=receivestruct.minsourcetranslationlevel)
    _aggregate!(sourcestruct, receivestruct.aggregationlist)
    _transfer!(sourcestruct, receivestruct)
    _disaggregate!(receivestruct)

    (sourcestruct.verbose || receivestruct.verbose) && println("---------------------------------")
    nothing
end
function _forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive, x)
    sourcestruct.xvector .= x
    _forward!(sourcestruct, receivestruct)
end

"""
    _transpose_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive, [y])

Store the result of the transposed matrix vector product in `sourcevector.xvector`

If no vector `y` is given, `receivestruct.bvector` is used as tight vector.
Otherwise, `receivestruct.bvector` is overwritten by the content of `y` before the operation is performed
"""
function _transpose_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)
    (sourcestruct.verbose || receivestruct.verbose) && @info "---------------------------------\n   Evaluate transpose forward operator\n---------------------------------------"

    _transpose_disaggregate!(receivestruct)
    _transpose_transfer!(sourcestruct, receivestruct)
    _transpose_aggregate!(sourcestruct, receivestruct.aggregationlist)       

    (sourcestruct.verbose || receivestruct.verbose) && println("-----------------------------------------")
    nothing
end
function _transpose_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive, y)
    receivestruct.bvector .= y
    _transpose_forward!(sourcestruct, receivestruct)
end

"""
    _adjoint_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)

Store the result of the adjoint matrix vector product in `sourcevector.xvector`

If no vector `y` is given, `receivestruct.bvector` is used as tight vector.
Otherwise, `receivestruct.bvector` is overwritten by the content of `y` before the operation is performed    
"""
function _adjoint_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive)
    (sourcestruct.verbose || receivestruct.verbose) && @info "---------------------------------\n   Evaluate adjoint forward operator\n---------------------------------------"

    _adjoint_disaggregate!(receivestruct)
    _adjoint_transfer!(sourcestruct, receivestruct)
    _adjoint_aggregate!(sourcestruct, receivestruct.aggregationlist)       

    (sourcestruct.verbose || receivestruct.verbose) && println("-----------------------------------------")
    nothing
end
function _adjoint_forward!(sourcestruct::MLFMMSource, receivestruct::MLFMMReceive, y)
    receivestruct.bvector .= y
    _adjoint_forward!(sourcestruct, receivestruct)
end


