
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
            reset = reset,
        )
        _muladd_or_mulreset!(
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
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(A.nodefarfields[parentnode]) .* A.phaseshifttoparent[sector, level+1]
        # A.nodefarfields[child] .=
        _muladd_or_mulreset!(
            A.nodefarfields[child],
            transpose_resamplemat,
            resamplemap.outputbuffer,
            reset = !A.nodeisfresh[child],
        )


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
        A.nodefarfields[child] .= _muladd_or_mulreset!(
            A.nodefarfields[child].buffer,
            adjoint_resamplemat,
            resamplemap.outputbuffer,
            reset = !A.nodeisfresh[child],
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
        for sector = 1:8
            # Threads.@threads for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end

        for node in aggregationlist[level]
            # Threads.@threads for node in aggregationlist[level]
            _transpose_aggregate_children!(A, node)
        end

        for sector = 1:8
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
function _aggregate_to_minlevel!(A::MLFMMSource, x; min_aggregationlevel::Integer = 0)
    A.buffer .= x
    _aggregate_to_minlevel!(A, min_aggregationlevel = min_aggregationlevel)
end
function _aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer = 0)
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
function _adjoint_aggregate_to_minlevel!(A::MLFMMSource; min_aggregationlevel::Integer = 0)
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
    min_aggregationlevel::Integer = 0,
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
