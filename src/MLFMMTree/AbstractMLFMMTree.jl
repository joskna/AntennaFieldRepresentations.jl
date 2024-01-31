abstract type AbstractMLFMMTree end

function (t::AbstractMLFMMTree)(index::Int)
    return error("Not implemented")
end

function ClusterTrees.DepthFirstIterator(tree::AbstractMLFMMTree, node::Int)
    return error("Not implemented")
end

function ClusterTrees.root(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function root(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function nearNodes(tree::AbstractMLFMMTree, node::Int, bufferboxes::Int)
    return error("Not implemented")
end

function nearnodeindices(tree::AbstractMLFMMTree, node::Int, bufferboxes::Int)
    return error("Not implemented")
end

function farNodes(tree::AbstractMLFMMTree, node::Int, bufferboxes::Int)
    return error("Not implemented")
end

function nodesatlevel(tree::AbstractMLFMMTree, level::Int)
    return error("Not implemented")
end

function samelevelnodes(tree::AbstractMLFMMTree, level::Int)
    return error("Not implemented")
end

function leaves(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function isleaf(tree::AbstractMLFMMTree, node::Int)
    return error("Not implemented")
end

function points(tree::AbstractMLFMMTree, node::Int)
    return error("Not implemented")
end

function level(tree::AbstractMLFMMTree, nodeid=ClusterTrees.root(tree)::Int)
    return error("Not implemented")
end

function levels(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function numberofnodes(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function center(tree::AbstractMLFMMTree, nodeid=ClusterTrees.root(tree)::Int)
    return error("Not implemented")
end

function halfsize(tree::AbstractMLFMMTree, nodeid=ClusterTrees.root(tree)::Int)
    return error("Not implemented")
end

function minhalfsize(tree::AbstractMLFMMTree)
    return error("Not implemented")
end

function ParentDepthFirstIterator(tree::AbstractMLFMMTree, node::Int)
    return error("Not implemented")
end

function transferingNodes(tree::AbstractMLFMMTree, node::Int, bufferboxes::Int)
    return error("Not implemented")
end

function nontransferingNodes(tree::AbstractMLFMMTree, node::Int, bufferboxes::Int)
    return error("Not implemented")
end
