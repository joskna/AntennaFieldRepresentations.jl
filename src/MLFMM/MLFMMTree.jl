using ClusterTrees: ClusterTrees

struct BoxData{N,T}
    values::Vector{Int}
    center::SVector{N,T}
    halfsize::T
    sector::Int
    level::Int
end

struct Node{B}
    data::B
    next_sibling::Int
    parent::Int
    first_child::Int
end

struct MLFMMTree{B,N,T} <: ClusterTrees.PointerBasedTrees.APBTree
    nodes::Vector{Node{B}}
    root::Int
    center::SVector{N,T}
    halfsize::T
    nodesatlevel::Vector{Vector{Int}}
end

(tree::MLFMMTree)(index::Integer) = tree.nodes[index]

function MLFMMTree(center::SVector{N,T}, halfsize::T) where {N,T}
    root = Node(BoxData(Int[], center, halfsize, 1, 1), 0, 0, 0)
    return MLFMMTree([root], 1, center, halfsize, [Int[]])
end

function MLFMMTree(
    center::SVector{N,T},
    points::AbstractArray{SVector{N,T},1},
    halfsize::T,
    minhalfsize::T,
) where {N,T}

    #ensure that halfsize is 2ᴺ ⋅ minhalfsize
    halfsize_powerof2 = log2(halfsize / minhalfsize)
    halfsize_powerof2 = maximum([0, Integer(ceil(halfsize_powerof2))])
    halfsize = 2^halfsize_powerof2 * minhalfsize
    tree = MLFMMTree(center, halfsize)
    addpoints!(tree, points, center, halfsize, minhalfsize)
    calculatenodesatlevel!(tree)

    for level ∈ tree.nodesatlevel
        sort!(level)
    end

    return tree
end

function calculatenodesatlevel!(tree::MLFMMTree)
    for node ∈ ClusterTrees.DepthFirstIterator(tree, root(tree))
        if level(tree, node) > length(nodesatlevel(tree))
            for i ∈ 1:(level(tree, node)-length(nodesatlevel(tree))-1)
                push!(tree.nodesatlevel, Int[])
            end
            push!(tree.nodesatlevel, [node])

        else
            append!(tree.nodesatlevel[MLFMMTrees.level(tree, node)], node)
        end
    end
end

function addpoints!(
    tree::MLFMMTree,
    points::AbstractArray{SVector{N,T},1},
    rootcenter::SVector{N,T},
    rootsize::T,
    smallestboxsize::T,
) where {N,T}
    for i ∈ eachindex(points)
        router = ClusterTrees.Octrees.Router(smallestboxsize, points[i])
        root_state = root(tree), rootcenter, rootsize, 1, 1
        ClusterTrees.update!(tree, root_state, i, router) do tree, node, data
            push!(data(tree, node).values, data)
        end
    end
end

function data(tree::MLFMMTree, node::Integer)
    return tree(node).data
end
function ClusterTrees.data(tree::MLFMMTree, node::Integer)
    return tree(node).data
end

function root(tree::MLFMMTree)
    return tree.root
end
function ClusterTrees.root(tree::MLFMMTree)
    return tree.root
end

function route!(tree::MLFMMTree, state, destination)
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
        child_sector = data(tree, child).sector
        child_pos = ClusterTrees.Octrees.hilbert_positions[sfc_state][child_sector+1] + 1
        child_level = data(tree, child).level
        target_pos < child_pos && break
        if child_sector == target_sector
            return child, target_center, target_size, target_sfc_state, child_level
        end
        pos = newpos
    end

    data = BoxData(Int[], target_center, target_size, target_sector, target_level)
    child = insert!(chds, data, pos)

    return child, target_center, target_size, target_sfc_state, target_level
end
function ClusterTrees.route!(tree::MLFMMTree, state, destination)
    return route!(tree::MLFMMTree, state, destination)
end

function Base.insert!(chd_itr::ClusterTrees.ChildIterator{MLFMMTree}, item, state)
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

function setfirstchild!(tree::MLFMMTree, nodeid, child)
    node = tree(nodeid)
    return tree.nodes[nodeid] = Node(node.data, node.next_sibling, node.parent, child)
end

function setnextsibling!(tree::MLFMMTree, nodeid, sibling)
    node = tree(nodeid)
    return tree.nodes[nodeid] = Node(node.data, sibling, node.parent, node.first_child)
end

function ClusterTrees.PointerBasedTrees.firstchild(tree::MLFMMTree, nodeid::Integer)
    return tree(nodeid).first_child
end

function parent(tree::MLFMMTree, node_idx::Integer)
    return tree(node_idx).parent
end
function ClusterTrees.parent(tree::MLFMMTree, node_idx::Integer)
    return tree(node_idx).parent
end

function ClusterTrees.PointerBasedTrees.nextsibling(tree::MLFMMTree, node_idx::Integer)
    return tree(node_idx).next_sibling
end

function center(tree::MLFMMTree, nodeid::Integer)
    return tree(nodeid).data.center
end

function halfsize(tree::MLFMMTree, nodeid::Integer)
    return tree(nodeid).data.halfsize
end

function level(tree::MLFMMTree, nodeid::Integer)
    return tree(nodeid).data.level
end


function isnear(
    center_a::AbstractVector,
    center_b::AbstractVector,
    halfsize::T,
    bufferboxes::Integer,
) where {T}
    Rvec = (center_a - center_b) ./ halfsize
    distancesquared = LinearAlgebra.dot(Rvec, Rvec)
    return distancesquared <= (bufferboxes + 1) * 12 * (1 + 10 * eps(T))
end

function isfar(
    center_a::AbstractVector,
    center_b::AbstractVector,
    halfsize::T,
    bufferboxes::Integer,
) where {T}
    return !isnear(center_a, center_b, halfsize, bufferboxes)
end

function samelevelnodes(tree::MLFMMTree, nodeid::Integer)
    level = level(tree, nodeid)
    return tree.nodesatlevel[level]
end

function getnontransferingNodes(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    return Iterators.filter(
        x -> (!istransferring(tree, node, x, bufferboxes)),
        samelevelnodes(tree, node),
    )
end

function gettransferingNodes(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    return Iterators.filter(
        x -> (istransferring(tree, node, x, bufferboxes)),
        samelevelnodes(tree, node),
    )
end

# function notransferMLFMM(tree::PBMLFMMTree, node::Integer, testnode::Integer, bufferboxes::Integer)
#     return !transferMLFMM(tree, node, testnode, bufferboxes)
# end

"""
    istransferring(sourcetree::MLFMMTree, [receivetree::MLFMMTree], sourcenode::Integer, receivenode::Integer, bufferboxes::Integer)

Returns true if an MLFMM transfer is planned between `tree(node)` and `tree(testnode)`.
"""
function istransferring(
    tree::MLFMMTree,
    node::Integer,
    testnode::Integer,
    bufferboxes::Integer,
)
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
function istransferring(
    sourcetree::MLFMMTree,
    receivetree::MLFMMTree,
    sourcenode::Integer,
    receivenode::Integer,
    bufferboxes::Integer,
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
        sourcenodeparent = parent(sourcetree, sourcenode)
        receivenodeparent = parent(receivetree, receivenode)

        isfar(
            center(sourcetree, sourcenodeparent),
            center(receivetree, receivenodeparent),
            halfsize(sourcetree, sourcenodeparent),
            bufferboxes,
        ) && return false
    end
    return true
end

function nearNodes(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    return Iterators.filter(
        x -> (isnear(
            center(tree, node),
            center(tree, x),
            halfsize(tree, node),
            bufferboxes,
        )),
        samelevelnodes(tree, node),
    )
end

function farNodes(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    return Iterators.filter(
        x ->
            (isfar(center(tree, node), center(tree, x), halfsize(tree, node), bufferboxes)),
        samelevelnodes(tree, node),
    )
end

function nearnodeindices(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    indices = Int[]
    for n ∈ nearNodes(tree, node, bufferboxes)
        append!(indices, points(tree, n))
    end
    return indices
end

function farnodeindices(tree::MLFMMTree, node::Integer, bufferboxes::Integer)
    indices = Int[]
    for n ∈ farNodes(tree, node, bufferboxes)
        append!(indices, points(tree, n))
    end
    return indices
end

function nodesatlevel(tree::MLFMMTree)
    return tree.nodesatlevel
end

function nodesatlevel(tree::MLFMMTree, level::Integer)
    return tree.nodesatlevel[level]
end

function points(tree::MLFMMTree, node::Integer)
    tree(node).first_child == 0 && return data(tree, node).values

    values = Int[]
    for i ∈ ClusterTrees.leaves(tree, node)
        values = [values; points(tree, i)]
    end
    return values
end

function leafs(tree::MLFMMTree, node::Integer = root(tree))
    return collect(Iterators.flatten(ClusterTrees.leaves(tree, node)))
end

function minhalfsize(tree::MLFMMTree)
    nodesatlowestlevel = nodesatlevel(tree)[end]
    nodesatlowestlevel == Int[] && return 0

    return data(tree, nodesatlowestlevel[begin]).halfsize
end


function levels(tree::MLFMMTree)
    return 1:length(nodesatlevel(tree))
end

function isleaf(tree::MLFMMTree, node::Integer)
    return tree(node).first_child == 0
end

function numberofnodes(tree::MLFMMTree)
    return length(tree.nodes)
end
struct ReverseDepthFirstIterator
    tree::MLFMMTree
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

function parents(tree::MLFMMTree, node::Integer)
    return collect(ParentUpwardsIterator(tree, node))
end

function findleafnode(tree::MLFMMTree, value::Int)
    for leaf ∈ leafs(tree)
        (value ∈ points(tree, leaf)) && return leaf
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

"""
    getboundingbox(points::Vector{SVector{D, F}})

Returns halfsize and center of bounding box of points. The halfsize is the half of the length of the edge of the bounding box.
"""
function getboundingbox(points::AbstractArray{SVector{D,F},1}) where {D,F}
    min_dim = Vector(points[1])
    max_dim = Vector(points[1])

    for i in eachindex(points)
        for j ∈ 1:D
            min_dim[j] = min_dim[j] < points[i][j] ? min_dim[j] : points[i][j]
            max_dim[j] = max_dim[j] > points[i][j] ? max_dim[j] : points[i][j]
        end
    end

    center = zeros(F, D)

    length_dim = zeros(F, D)
    for j ∈ 1:D
        length_dim[j] = max_dim[j] - min_dim[j]
        center[j] = (max_dim[j] + min_dim[j]) / F(2.0)
    end

    halflength = maximum(length_dim) / F(2.0)

    return SVector{D}(center), halflength
end
