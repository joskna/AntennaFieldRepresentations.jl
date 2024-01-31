struct MLFMMTree{N,D,T} <: MLFMMTrees.AbstractMLFMMTree
    tree::MLFMMTrees.PBMLFMMTree{N,D,T}
end

function MLFMMTree(
    center::SVector{D,T}, points::Vector{SVector{D,T}}, halfsize::T, minhalfsize::T
) where {D,T}
    return MLFMMTree(MLFMMTrees.PBMLFMMTree(center, points, halfsize, minhalfsize))
end

function tree(tree::MLFMMTree)
    return tree.tree
end

function (t::MLFMMTree)(index::Int)
    return MLFMMTrees.tree(t)(index)
end

function ClusterTrees.DepthFirstIterator(tree::MLFMMTree, node::Int)
    return ClusterTrees.DepthFirstIterator(MLFMMTrees.tree(tree), node)
end

function ClusterTrees.root(tree::MLFMMTree)
    return ClusterTrees.root(MLFMMTrees.tree(tree))
end

function root(tree::MLFMMTree)
    return MLFMMTrees.root(MLFMMTrees.tree(tree))
end

function nontransferingNodes(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.nontransferingNodes(MLFMMTrees.tree(tree), node, bufferboxes)
end

function transferingNodes(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.transferingNodes(MLFMMTrees.tree(tree), node, bufferboxes)
end

function samelevelnodes(tree::MLFMMTree, level::Int)
    return samelevelnodes(MLFMMTrees.tree(tree), level)
end

function leafs(tree::MLFMMTree)
    return MLFMMTrees.leafs(MLFMMTrees.tree(tree))
end

function level(tree::MLFMMTree, nodeid::Int)
    return MLFMMTrees.level(MLFMMTrees.tree(tree), nodeid)
end

function center(tree::MLFMMTree, nodeid=ClusterTrees.root(tree)::Int)
    return MLFMMTrees.center(MLFMMTrees.tree(tree), nodeid)
end

function halfsize(tree::MLFMMTree, nodeid=ClusterTrees.root(tree)::Int)
    return MLFMMTrees.halfsize(MLFMMTrees.tree(tree), nodeid)
end

function points(tree::MLFMMTree, node::Int)
    return MLFMMTrees.points(MLFMMTrees.tree(tree), node)
end

function minhalfsize(tree::MLFMMTree)
    return MLFMMTrees.minhalfsize(MLFMMTrees.tree(tree))
end

function nearnodeindices(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.nearnodeindices(MLFMMTrees.tree(tree), node, bufferboxes)
end

function farnodeindices(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.farnodeindices(MLFMMTrees.tree(tree), node, bufferboxes)
end

function nodesatlevel(tree::MLFMMTree, level::Int)
    return MLFMMTrees.nodesatlevel(MLFMMTrees.tree(tree), level)
end

function levels(tree::MLFMMTree)
    return MLFMMTrees.levels(MLFMMTrees.tree(tree))
end

function isleaf(tree::MLFMMTree, node::Int)
    return MLFMMTrees.isleaf(MLFMMTrees.tree(tree), node)
end

function numberofnodes(tree::MLFMMTree)
    return MLFMMTrees.numberofnodes(MLFMMTrees.tree(tree))
end

function ClusterTrees.children(tree::MLFMMTree, node::Int)
    return ClusterTrees.children(MLFMMTrees.tree(tree), node)
end

function ClusterTrees.parent(tree::MLFMMTree, node_idx::Int)
    return ClusterTrees.parent(MLFMMTrees.tree(tree), node_idx)
end

function parents(tree::MLFMMTree, node_idx::Int)
    return MLFMMTrees.parents(MLFMMTrees.tree(tree), node_idx)
end

function nearNodes(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.nearNodes(MLFMMTrees.tree(tree), node, bufferboxes)
end

function farNodes(tree::MLFMMTree, node::Int, bufferboxes::Int)
    return MLFMMTrees.farNodes(MLFMMTrees.tree(tree), node, bufferboxes)
end

function findleafnode(tree::MLFMMTree, value::Int)
    return MLFMMTrees.findleafnode(MLFMMTrees.tree(tree), value)
end

function ParentDepthFirstIterator(tree::MLFMMTree, node::Int)
    return MLFMMTrees.ParentDepthFirstIterator(MLFMMTrees.tree(tree), node)
end

import StaticArrays: MVector

"""
    getboundingbox(points::Vector{SVector{D, F}})

Returns halfsize and center of bounding box of points. The halfsize is the half of the length of the edge of the bounding box.
"""
function getboundingbox(points::AbstractArray{SVector{D,F},1}) where {D,F}
    min_dim = Vector(points[1])
    max_dim = Vector(points[1])

    for i ∈ 1:length(points)
        for j ∈ 1:D
            min_dim[j] = min_dim[j] < points[i][j] ? min_dim[j] : points[i][j]
            max_dim[j] = max_dim[j] > points[i][j] ? max_dim[j] : points[i][j]
        end
    end

    center = MVector{D}(zeros(F, D))

    length_dim = zeros(F, D)
    for j ∈ 1:D
        length_dim[j] = max_dim[j] - min_dim[j]
        center[j] = (max_dim[j] + min_dim[j]) / F(2.0)
    end

    halflength = maximum(length_dim) / F(2.0)

    return center, halflength
end