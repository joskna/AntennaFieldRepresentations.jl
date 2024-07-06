using StaticArrays
struct Data{D,T}
    values::Vector{Int}
    center::SVector{D,T}
    halfsize::T
    sector::Int
    level::Int
end

# in comparison to ClusterTrees in https://github.com/krcools/ClusterTrees.jl/blob/df3b0ecb7531097a9a7d8bcc4c009f4ebd0078f3/src/pointerbasedtrees.jl#L9 we removed 
# the num_children field as it wasn't being set
struct Node{N}
    data::N
    next_sibling::Int
    parent::Int
    first_child::Int
end

struct PBMLFMMTree{N,D,T} <: ClusterTrees.PointerBasedTrees.APBTree
    nodes::Vector{Node{N}}
    root::Int
    center::SVector{D,T}
    halfsize::T
    nodesatlevel::Vector{Vector{Int}}
end

(tree::PBMLFMMTree)(index::Int) = tree.nodes[index]

function PBMLFMMTree(center::SVector{D,T}, halfsize::T) where {D,T}
    root = Node(Data(Int[], center, halfsize, 1, 1), 0, 0, 0)
    return PBMLFMMTree([root], 1, center, halfsize, [Int[]])
end

function PBMLFMMTree(
    center::SVector{D,T},
    points::AbstractArray{SVector{D,T},1},
    halfsize::T,
    minhalfsize::T,
) where {D,T}

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2 = log2(halfsize / minhalfsize)
    halfsize_powerof2 = maximum([0, Integer(ceil(halfsize_powerof2))])
    halfsize = 2^halfsize_powerof2 * minhalfsize
    tree = PBMLFMMTree(center, halfsize)
    addpoints!(tree, points, center, halfsize, minhalfsize)
    calculatenodesatlevel!(tree)

    for level ∈ tree.nodesatlevel
        sort!(level)
    end

    return tree
end

function calculatenodesatlevel!(tree::PBMLFMMTree)
    for node ∈ ClusterTrees.DepthFirstIterator(tree, MLFMMTrees.root(tree))
        if MLFMMTrees.level(tree, node) > length(nodesatlevel(tree))
            for i ∈ 1:(MLFMMTrees.level(tree, node)-length(nodesatlevel(tree))-1)
                push!(tree.nodesatlevel, Int[])
            end
            push!(tree.nodesatlevel, [node])

        else
            append!(tree.nodesatlevel[MLFMMTrees.level(tree, node)], node)
        end
    end
end
function addpoints!(
    tree::PBMLFMMTree,
    points::AbstractArray{SVector{D,T},1},
    rootcenter::SVector{D,T},
    rootsize::T,
    smallestboxsize::T,
) where {D,T}
    for i ∈ eachindex(points)
        router = ClusterTrees.Octrees.Router(smallestboxsize, points[i])
        root_state = MLFMMTrees.root(tree), rootcenter, rootsize, 1, 1
        ClusterTrees.update!(tree, root_state, i, router) do tree, node, data
            push!(ClusterTrees.data(tree, node).values, data)
        end
    end
end

function ClusterTrees.data(tree::PBMLFMMTree, node::Int)
    return tree(node).data
end

function root(tree::PBMLFMMTree)
    return tree.root
end

function ClusterTrees.root(tree::PBMLFMMTree)
    return MLFMMTrees.root(tree)
end

function ClusterTrees.route!(tree::PBMLFMMTree, state, destination)
    point = destination.target_point
    smallest_box_size = destination.smallest_box_size

    nodeid, center, size, sfc_state, level = state
    size <= smallest_box_size && return state
    target_sector, target_center, target_size =
        ClusterTrees.Octrees.sector_center_size(point, center, size)
    target_pos = ClusterTrees.Octrees.hilbert_positions[sfc_state][target_sector+1] + 1
    target_sfc_state = ClusterTrees.Octrees.hilbert_states[sfc_state][target_sector+1] + 1
    target_level = level + 1

    chds = ClusterTrees.children(tree, nodeid)
    pos = ClusterTrees.start(chds)
    while !ClusterTrees.done(chds, pos)
        child, newpos = ClusterTrees.next(chds, pos)
        child_sector = ClusterTrees.data(tree, child).sector
        child_pos = ClusterTrees.Octrees.hilbert_positions[sfc_state][child_sector+1] + 1
        child_level = ClusterTrees.data(tree, child).level
        target_pos < child_pos && break
        if child_sector == target_sector
            return child, target_center, target_size, target_sfc_state, child_level
        end
        pos = newpos
    end

    data = Data(Int[], target_center, target_size, target_sector, target_level)
    child = insert!(chds, data, pos)

    return child, target_center, target_size, target_sfc_state, target_level
end

function Base.insert!(chd_itr::ClusterTrees.ChildIterator{<:PBMLFMMTree}, item, state)
    prev, next = state
    parent = chd_itr.node

    tree = chd_itr.tree
    push!(tree.nodes, Node(item, next, parent, 0))
    this = lastindex(tree.nodes)
    if prev < 1
        setfirstchild!(tree, parent, this)
    else
        setnextsibling!(tree, prev, this)
    end
    return this
end

function setfirstchild!(tree::PBMLFMMTree, nodeid, child)
    node = tree(nodeid)
    return tree.nodes[nodeid] = Node(node.data, node.next_sibling, node.parent, child)
end

function setnextsibling!(tree::PBMLFMMTree, nodeid, sibling)
    node = tree(nodeid)

    return tree.nodes[nodeid] = Node(node.data, sibling, node.parent, node.first_child)
end

function ClusterTrees.PointerBasedTrees.firstchild(tree::PBMLFMMTree, nodeid::Int)
    return tree(nodeid).first_child
end

function ClusterTrees.parent(tree::PBMLFMMTree, node_idx::Int)
    return tree(node_idx).parent
end

function ClusterTrees.PointerBasedTrees.nextsibling(tree::PBMLFMMTree, node_idx::Int)
    return tree(node_idx).next_sibling
end

function center(tree::PBMLFMMTree, nodeid::Int)
    return tree(nodeid).data.center
end

function halfsize(tree::PBMLFMMTree, nodeid::Int)
    return tree(nodeid).data.halfsize
end

function level(tree::PBMLFMMTree, nodeid::Int)
    return tree(nodeid).data.level
end

function isnear(
    center_a::SVector{D,T},
    center_b::SVector{D,T},
    halfsize::T,
    bufferboxes::Int,
) where {D,T}
    Rvec = (center_a - center_b) ./ halfsize
    distancesquared = LinearAlgebra.dot(Rvec, Rvec)
    return distancesquared <= (bufferboxes + 1) * 12 * (1 + 10 * eps(T))
end

function isfar(
    center_a::SVector{D,T},
    center_b::SVector{D,T},
    halfsize::T,
    bufferboxes::Int,
) where {D,T}
    return !isnear(center_a, center_b, halfsize, bufferboxes)
end

function samelevelnodes(tree::PBMLFMMTree, nodeid::Int)
    level = MLFMMTrees.level(tree, nodeid)
    return tree.nodesatlevel[level]
end

function nontransferingNodes(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    return Iterators.filter(
        x -> (notransferMLFMM(tree, node, x, bufferboxes)),
        MLFMMTrees.samelevelnodes(tree, node),
    )
end

function transferingNodes(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    return Iterators.filter(
        x -> (transferMLFMM(tree, node, x, bufferboxes)),
        MLFMMTrees.samelevelnodes(tree, node),
    )
end

function notransferMLFMM(tree::PBMLFMMTree, node::Int, testnode::Int, bufferboxes::Int)
    return !transferMLFMM(tree, node, testnode, bufferboxes)
end

function transferMLFMM(tree::PBMLFMMTree, node::Int, testnode::Int, bufferboxes::Int)
    MLFMMTrees.level(tree, node) < 3 && return false

    if isnear(center(tree, node), center(tree, testnode), halfsize(tree, node), bufferboxes)
        return false
    elseif MLFMMTrees.level(tree, node) > 3
        nodeparent = ClusterTrees.parent(tree, node)
        testnodeparent = ClusterTrees.parent(tree, testnode)

        isfar(
            center(tree, nodeparent),
            center(tree, testnodeparent),
            halfsize(tree, nodeparent),
            bufferboxes,
        ) && return false
    end
    return true
end

#TODO: better name for this function
function transferMLFMM(
    sourcetree::PBMLFMMTree,
    receivetree::PBMLFMMTree,
    sourcenode::Int,
    receivenode::Int,
    bufferboxes::Int,
)

    if isnear(
        center(sourcetree, sourcenode),
        center(receivetree, receivenode),
        halfsize(sourcetree, sourcenode),
        bufferboxes,
    )
        if isleaf(receivetree, receivenode)
            ErrorException("Probes too close to source tree! Cannot plan transfers.")
        else
            return false
        end
    elseif MLFMMTrees.level(sourcetree, sourcenode) > 1 &&
           MLFMMTrees.level(receivetree, receivenode) > 1
        sourcenodeparent = ClusterTrees.parent(sourcetree, sourcenode)
        receivenodeparent = ClusterTrees.parent(receivetree, receivenode)

        isfar(
            center(sourcetree, sourcenodeparent),
            center(receivetree, receivenodeparent),
            halfsize(sourcetree, sourcenodeparent),
            bufferboxes,
        ) && return false
    end
    return true
end

function nearNodes(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    return Iterators.filter(
        x -> (isnear(
            center(tree, node),
            center(tree, x),
            halfsize(tree, node),
            bufferboxes,
        )),
        MLFMMTrees.samelevelnodes(tree, node),
    )
end

function farNodes(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    return Iterators.filter(
        x ->
            (isfar(center(tree, node), center(tree, x), halfsize(tree, node), bufferboxes)),
        MLFMMTrees.samelevelnodes(tree, node),
    )
end

function nearnodeindices(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    indices = Int[]
    for n ∈ MLFMMTrees.nearNodes(tree, node, bufferboxes)
        append!(indices, MLFMMTrees.points(tree, n))
    end
    return indices
end

function farnodeindices(tree::PBMLFMMTree, node::Int, bufferboxes::Int)
    indices = Int[]
    for n ∈ MLFMMTrees.farNodes(tree, node, bufferboxes)
        append!(indices, MLFMMTrees.points(tree, n))
    end
    return indices
end

function nodesatlevel(tree::PBMLFMMTree)
    return tree.nodesatlevel
end

function nodesatlevel(tree::PBMLFMMTree, level::Int)
    return tree.nodesatlevel[level]
end

function points(tree::PBMLFMMTree, node::Int)
    tree(node).first_child == 0 && return ClusterTrees.data(tree, node).values

    values = Int[]
    for i ∈ ClusterTrees.leaves(tree, node)
        values = [values; MLFMMTrees.points(tree, i)]
    end
    return values
end

function leafs(tree::PBMLFMMTree, node::Int = ClusterTrees.root(tree))
    return collect(Iterators.flatten(ClusterTrees.leaves(tree, node)))
end

function minhalfsize(tree::PBMLFMMTree)
    nodesatlowestlevel = MLFMMTrees.nodesatlevel(tree)[end]
    nodesatlowestlevel == Int[] && return 0

    return ClusterTrees.data(tree, nodesatlowestlevel[begin]).halfsize
end

function levels(tree::PBMLFMMTree)
    return 1:length(MLFMMTrees.nodesatlevel(tree))
end

function isleaf(tree::PBMLFMMTree, node::Int)
    return tree(node).first_child == 0
end

function numberofnodes(tree::PBMLFMMTree)
    return length(tree.nodes)
end
struct ReverseDepthFirstIterator
    tree::PBMLFMMTree
    node::Int
end

Base.IteratorSize(::ReverseDepthFirstIterator) = Base.SizeUnknown()

function Base.iterate(itr::ReverseDepthFirstIterator)
    # chitr = ClusterTrees.children(itr.tree, itr.node)
    # stack::Vector{Int} = collect(chitr)

    stack::Vector{Int} = collect(ClusterTrees.children(itr.tree, itr.node))

    return itr.node, stack
end

function Base.iterate(itr::ReverseDepthFirstIterator, stack::Vector{Int})
    isempty(stack) && return nothing

    # node = stack[begin]
    # popfirst!(stack)
    node = popfirst!(stack)
    prepend!(stack, collect(ClusterTrees.children(itr.tree, node)))

    return node, stack
end

struct ParentUpwardsIterator{T}
    tree::T
    node::Int
end
Base.IteratorSize(::ParentUpwardsIterator) = Base.SizeUnknown()

function Base.iterate(itr::ParentUpwardsIterator)
    parent = ClusterTrees.parent(itr.tree, itr.node)
    if parent == 0
        return 0, Int[]
    end
    stack = Int[parent]
    return parent, stack
end

function Base.iterate(itr::ParentUpwardsIterator, stack)
    isempty(stack) && return nothing

    node = stack[begin]
    popfirst!(stack)
    parent = ClusterTrees.parent(itr.tree, node)
    if parent == 0
        return nothing
    else
        pushfirst!(stack, parent)
    end

    return parent, stack
end

function parents(tree::PBMLFMMTree, node::Int)
    return collect(ParentUpwardsIterator(tree, node))
end
# function level(tree::PBMLFMMTree, node::Int)
#     levels=MLFMMTrees.levels(tree)
#     for level in levels
#         a= findfirst( x-> x==node, MLFMMTrees.nodesatlevel(tree, level))
#         if a !==Nothing
#             return level
#         end
#     end
#     return 0
# end

function findleafnode(tree::PBMLFMMTree, value::Int)
    for leaf ∈ MLFMMTrees.leafs(tree)
        (value ∈ MLFMMTrees.points(tree, leaf)) && return leaf
    end
    return 0
end






struct NodeInformation
    info::Union{Nothing,Tuple{Int,Tuple{Int,Int}}}
end

function node(next::NodeInformation)
    return next.info[1]
end

function state(next::NodeInformation)
    return next.info[2]
end

function Base.isnothing(x::NodeInformation)
    return isnothing(x.info)
end

struct StackElement{T,N}
    chitr::ClusterTrees.ChildIterator{T,N}
    info::NodeInformation
end

function StackElement(chitr::ClusterTrees.ChildIterator, info::Tuple{Int,Tuple{Int,Int}})
    return StackElement(chitr, NodeInformation(info))
end

function StackElement(chitr::ClusterTrees.ChildIterator, info::Nothing)
    return StackElement(chitr, NodeInformation(info))
end

function childreniterator(s::StackElement{T,N}) where {T,N}
    return s.chitr
end

function information(s::StackElement{T,N}) where {T,N}
    return s.info
end

struct DepthFirstIterator{T,N}
    tree::T
    node::N
end

function Base.iterate(itr::DepthFirstIterator{T,N}) where {T,N}
    chitr = ClusterTrees.children(itr.tree, itr.node)
    stack = StackElement{T,N}[StackElement(chitr, iterate(chitr))]
    return iterate(itr, stack)
end

function Base.iterate(itr::DepthFirstIterator{T,N}, stack) where {T,N}
    isempty(stack) && return nothing
    while true
        info = information(last(stack))
        if !isnothing(info)
            n = node(info)
            chitr = ClusterTrees.children(itr.tree, n)
            push!(stack, StackElement(chitr, iterate(chitr)))
        else
            pop!(stack)
            isempty(stack) && return itr.node, stack

            chitr = childreniterator(last(stack))
            info = information(last(stack))
            n = node(info)
            s = state(info)
            stack[end] = StackElement(chitr, iterate(chitr, s))

            return n, stack
        end
    end
end
