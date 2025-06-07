

"""
    _disaggregate_leafnodes!(receivestruct::MLFMMReceive)

Test leaf node patterns with all test functions of corresponding leafnode and store result in `receivestruct.bvector`
"""
function _disaggregate_leafnodes!(receivestruct::MLFMMReceive)
    receivetree = receivestruct.tree

    # Threads.@threads for leafnode in receivestruct.leafnodeindices
    for leafnode in receivestruct.leafnodeindices
        !(receivestruct.nodeisoccupied[leafnode]) && continue

        pws = receivestruct.nodespectra[leafnode]

        for probeindex::Int in receivetree(leafnode).data.values::Vector{Int}
            ff = receivestruct.testfunctionfarfields[probeindex]

            receivestruct.buffer[probeindex] = udot(pws, ff)
            # receivestruct.buffer[probeindex] = sum(pws .* ff)
        end
    end

end


"""
    _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)

Perform ajoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate_leafnodes!`
"""
function _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)
    # function _transpose_disaggregate_leafnodes!(A::MLFMMReceive)
    receivetree = A.tree

    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        reset = true
        for functionindex in receivetree(leafnode).data.values
            if !reset
                _eθ(A.nodespectra[leafnode]) .+=
                    A.buffer[functionindex] .*
                    conj.(_eθ(A.testfunctionfarfields[functionindex]))
                _eϕ(A.nodespectra[leafnode]) .+=
                    A.buffer[functionindex] .*
                    conj.(_eϕ(A.testfunctionfarfields[functionindex]))
            else
                _eθ(A.nodespectra[leafnode]) .=
                    A.buffer[functionindex] .*
                    conj.(_eθ(A.testfunctionfarfields[functionindex]))
                _eϕ(A.nodespectra[leafnode]) .=
                    A.buffer[functionindex] .*
                    conj.(_eϕ(A.testfunctionfarfields[functionindex]))
                reset = false
            end
        end
    end

end

"""
    _transpose_disaggregate_leafnodes!(A::MLFMMReceive)

Perform transposed operation to `disaggregate_leafnodes!`
"""
function _transpose_disaggregate_leafnodes!(A::MLFMMReceive)
    # function _adjoint_disaggregate_leafnodes!(A::MLFMMReceive)
    receivetree = A.tree

    # Threads.@threads for leafnode in A.leafnodeindices
    for leafnode in A.leafnodeindices
        !(A.nodeisoccupied[leafnode]) && continue
        reset = true
        for functionindex::Int in receivetree(leafnode).data.values::Vector{Int}
            A.nodespectra[leafnode].buffer .= _muladd_or_mulreset!(
                A.nodespectra[leafnode],
                A.testfunctionfarfields[functionindex],
                A.buffer[functionindex],
                reset = reset,
            )
            reset = false
        end
    end

end

"""
    _disaggregate_children!(receivestruct, parentnode)

Disaggregate pattern of parentnode to all its children and store the resulting patterns in `receivestruct.nodespectra[child]`
"""
function _disaggregate_children!(receivestruct::MLFMMReceive, parentnode)

    receivetree = receivestruct.tree


    lvl = AntennaFieldRepresentations.level(receivetree, parentnode)

    resamplemap = receivestruct.levelresamplemaps[lvl]
    adjoint_resamplemat = adjoint(resamplemap)
    # reset = true

    for child in children(receivetree, parentnode)
        !(receivestruct.nodeisoccupied[child]) && continue

        sector = receivetree.nodes[child].data.sector + 1

        view(resamplemap.outputbuffermat, :, :, 1) .=
            _eθ(receivestruct.nodespectra[parentnode]) .*
            (receivestruct.phaseshifttoparent[sector, lvl+1])
        view(resamplemap.outputbuffermat, :, :, 2) .=
            _eϕ(receivestruct.nodespectra[parentnode]) .*
            (receivestruct.phaseshifttoparent[sector, lvl+1])
        # receivestruct.nodespectra[child] .=
        # _adjoint_resamplematrix!(view(resamplemap.outputbuffermat, :, :, 1), _eθ(receivestruct.nodespectra[child]), resamplemap; reset=!receivestruct.nodeisfresh[child])
        # _adjoint_resamplematrix!(view(resamplemap.outputbuffermat, :, :, 2), _eϕ(receivestruct.nodespectra[child]), resamplemap; reset=!receivestruct.nodeisfresh[child])

        _muladd_or_mulreset!(
            receivestruct.nodespectra[child],
            adjoint_resamplemat,
            resamplemap.outputbuffer,
            reset = !receivestruct.nodeisfresh[child],
        )
        receivestruct.nodeisfresh[child] = true
    end
end

"""
    _transpose_disaggregate_children!(A, parentnode)

Perform transposed operation to `disaggrgate_children!`
"""
function _transpose_disaggregate_children!(A, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)

    resamplemap = A.levelresamplemaps[level]
    reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        mul!(resamplemap.outputbuffer, resamplemap, A.nodespectra[child])

        _muladd_or_mulreset!(
            _eθ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            (A.phaseshifttoparent[sector, level+1]),
            reset = reset,
        )
        _muladd_or_mulreset!(
            _eϕ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            (A.phaseshifttoparent[sector, level+1]),
            reset = reset,
        )

        reset = false
    end
end

"""
    _adjoint_disaggregate_children!(A, parentnode)

Perform adjoint operation to `disaggrgate_children!`
"""
function _adjoint_disaggregate_children!(A, parentnode)
    tree = A.tree

    level = AntennaFieldRepresentations.level(tree, parentnode)


    # resamplemap = conj.(A.levelresamplemaps[level]) -> since resamplemap is real, conj.(conj.(A.levelresamplemaps[level]) == A.levelresamplemaps[level])
    resamplemap = A.levelresamplemaps[level]
    reset = true

    for child in children(tree, parentnode)
        !(A.nodeisoccupied[child]) && continue

        sector = tree.nodes[child].data.sector + 1

        mul!(resamplemap.outputbuffer, resamplemap, A.nodespectra[child])

        _muladd_or_mulreset!(
            _eθ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 1),
            (A.phaseshifttoparent[sector, level+1]),
            reset = reset,
        )
        _muladd_or_mulreset!(
            _eϕ(A.nodespectra[parentnode]),
            view(resamplemap.outputbuffermat, :, :, 2),
            (A.phaseshifttoparent[sector, level+1]),
            reset = reset,
        )

        reset = false
    end
end

"""
    _disaggregate_to_disaggregationslist!(receivestruct::MLFMMReceive)

Perform all disaggregations according to `receivestruct.disaggregationlist`.
The `receivestruct.disaggregationlist` stores all nodes at each level which shall perform a disaggregation.
"""
function _disaggregate_to_disaggregationslist!(receivestruct::MLFMMReceive)
    receivestruct.verbose && @info "Disaggregate node spectra "


    for level in eachindex(receivestruct.disaggregationlist)
        # for level in levels(receivestruct.tree)
        receivestruct.disaggregationlist[level] == [] && continue
        for sector = 1:8
            # Threads.@threads for sector = 1:8
            conj!(receivestruct.phaseshifttoparent[sector, level+1])
        end

        for node in receivestruct.disaggregationlist[level]
            # Threads.@threads for node in receivestruct.disaggregationlist[level]
            !(receivestruct.nodeisoccupied[node]) && continue
            _disaggregate_children!(receivestruct, node)
        end

        for sector = 1:8
            # Threads.@threads for sector = 1:8
            conj!(receivestruct.phaseshifttoparent[sector, level+1])
        end

    end

    _disaggregate_leafnodes!(receivestruct)

end

"""
    _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive, [y::AbstractVector])

Perform adjoint operation (i.e., complex conjugate of transposed operation) of `_disaggregate!`#

If no vector `y` is given, the content in `A.buffer` is used for adjoint disaggregation. 
Otherwise, `A.buffer` is overwritten by `y` before adjoint disaggregation.
"""
function _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive, y::AbstractVector)
    A.buffer .= y
    _adjoint_disaggregate_to_disaggregationslist!(A)
end
function _adjoint_disaggregate_to_disaggregationslist!(A::MLFMMReceive)
    A.verbose && @info "Adjoint disaggregate node spectra "

    _adjoint_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))

        for node in A.disaggregationlist[level]
            # Threads.@threads for node in A.disaggregationlist[level]
            !(A.nodeisoccupied[node]) && continue
            _adjoint_disaggregate_children!(A, node)
            # _transpose_disaggregate_children!(A, node)
        end

    end
    # _adjoint_disaggregate_leafnodes!(A)

end

"""
    _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive, [y::AbstractVector])

    Perform transposed operation of `_disaggregate!`

If no vector `y` is given, the content in `A.buffer` is used for transposed disaggregation. 
Otherwise, `A.buffer` is overwritten by `y` before transposed disaggregation.
"""
function _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive, y::AbstractVector)
    A.buffer .= y
    _transpose_disaggregate_to_disaggregationslist!(A)
end
function _transpose_disaggregate_to_disaggregationslist!(A::MLFMMReceive)
    A.verbose && @info "Transpose disaggregate node spectra "

    _transpose_disaggregate_leafnodes!(A)

    for level in reverse(eachindex(A.disaggregationlist))
        A.disaggregationlist[level] == [] && continue

        # Threads.@threads for sector = 1:8
        for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end


        # Threads.@threads for node in A.disaggregationlist[level]
        for node in A.disaggregationlist[level]
            !(A.nodeisoccupied[node]) && continue
            _transpose_disaggregate_children!(A, node)
        end

        # Threads.@threads for sector = 1:8
        for sector = 1:8
            conj!(A.phaseshifttoparent[sector, level+1])
        end

    end

end
