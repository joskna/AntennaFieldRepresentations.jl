##### 
#= The methods in this file are derived from the package BEAST.jl ( https://github.com/krcools/BEAST.jl ).
BEAST.jl is licenced under the MIT "Expat" Licence with the following licence text:

    Copyright (c) 2017: Kristof Cools.

    Permission is hereby granted, free of charge, to any person obtaining a copy

    of this software and associated documentation files (the "Software"), to deal

    in the Software without restriction, including without limitation the rights

    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell

    copies of the Software, and to permit persons to whom the Software is

    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all

    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR

    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,

    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE

    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER

    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,

    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

    SOFTWARE.
=#

struct typeindicator{T} end

function quaddata(op::BEAST.MWFarField, rs, els, qs::BEAST.SingleNumQStrat)
    return quadpoints(rs, els, (qs.quad_rule,))
end

function quadrule(
    op::BEAST.MWFarField3D,
    refspace,
    p,
    y,
    q,
    el,
    qdata,
    qs::BEAST.SingleNumQStrat,
)
    return qdata[1, q]
end
function quadrule(
    op::BEAST.MWDoubleLayerFarField3D,
    refspace,
    p,
    y,
    q,
    el,
    qdata,
    qs::BEAST.SingleNumQStrat,
)
    return qdata[1, q]
end

function localpotential(
    op,
    points,
    basis;
    t::typeindicator{T} = typeindicator{SVector{3,ComplexF64}}(),
    quadstrat = BEAST.SingleNumQStrat(7),
    verbose = false,
) where {T}
    ff = zeros(T, length(points), BEAST.numfunctions(basis))
    # ff = Vector{Array{type}}(undef,numfunctions(basis))
    # fill!(ff,fill!(Array{type}(undef,size(points)), type([0;0;0])))
    # store(v, p, n) = (ff[n][p] += v)
    store(v, m, n) = (ff[m, n] += v)
    localpotential!(
        store,
        op,
        points,
        basis;
        t = t,
        quadstrat = quadstrat,
        verbose = verbose,
    )
    return ff
end

# get the far-field for every ansatz function separately as a matrix
# far-fields are sampled at same points

function localpotential!(
    store,
    op,
    points,
    basis;
    t::typeindicator{T} = typeindicator{SVector{3,ComplexF64}}(),
    quadstrat = BEAST.SingleNumQStrat(7),
    verbose = false,
) where {T}
    verbose && @info "Multi-threaded assembly: ($(Threads.nthreads()) threads)"
    # z = zeros(type, length(points))

    els, ad = assemblydata(basis)
    rs = refspace(basis)
    Nthreads = Threads.nthreads()
    # zlocal = Array{T}(undef, BEAST.numfunctions(rs))
    zlocal = [Array{T}(undef, BEAST.numfunctions(rs)) for _ = 1:Nthreads]
    qdata = quaddata(op, rs, els, quadstrat)

    # pointsiterator = verbose ? ProgressBar(enumerate(points)) : enumerate(points)

    # pointsiterator = collect(enumerate(points))
    # for (p, y) ∈ pointsiterator
    Threads.@threads for p ∈ eachindex(points)
        threadid = Threads.threadid()
        y = points[p]
        for (q, el) ∈ enumerate(els)
            fill!(zlocal[threadid], zero(T))
            qr = quadrule(op, rs, p, y, q, el, qdata, quadstrat)
            BEAST.farfieldlocal!(zlocal[threadid], op, rs, y, el, qr)

            # assemble from local contributions
            for (r, z) ∈ enumerate(zlocal[threadid])
                for (n, b) ∈ ad[q, r]
                    store(z * b, p, n)
                end
            end
        end
    end
end
function individualcartesianfarfields(
    basisfunctions::SurfaceCurrentDensity,
    pts;
    # quadstrat = BEAST.SingleNumQStrat(1),
    quadstrat = BEAST.SingleNumQStrat(4),
    verbose = false,
)
    return individualcartesianfarfields(
        basisfunctions,
        pts,
        getwavenumber(basisfunctions);
        quadstrat = quadstrat,
        verbose = verbose,
    )
end
function individualcartesianfarfields(
    basisfunctions::SurfaceCurrentDensity{Radiated,Electric,S,C},
    pts,
    k0::T;
    # quadstrat = BEAST.SingleNumQStrat(1),
    quadstrat = BEAST.SingleNumQStrat(4),
    verbose = false,
) where {T<:Real,S<:BEAST.Space{T},C}
    funspace = functionspace(basisfunctions)
    factor = complex(zero(T), -k0) * Z₀ / (4π)
    ffs = localpotential(
        MWFarField3D(; wavenumber = k0),
        pts,
        funspace;
        t = typeindicator{SVector{3,Complex{T}}}(),
        quadstrat = quadstrat,
        verbose = verbose,
    )
    ffs .*= factor
    return ffs
end
function individualcartesianfarfields(
    basisfunctions::SurfaceCurrentDensity{Radiated,Magnetic,S,C},
    pts,
    k0::T;
    # quadstrat = BEAST.SingleNumQStrat(1),
    quadstrat = BEAST.SingleNumQStrat(4),
    verbose = false,
) where {T<:Real,S<:BEAST.Space{T},C}
    funspace = functionspace(basisfunctions)
    #TODO: check this factor
    # factor = complex(zero(R), -k0) / (4π)
    factor = complex(zero(T), -k0) * Z₀ / (4π) # same factor as for electric currents
    # factor=1
    ffs = localpotential(
        BEAST.MWDoubleLayerFarField3D(; wavenumber = k0),
        pts,
        funspace;
        t = typeindicator{SVector{3,Complex{R}}}(),
        quadstrat = quadstrat,
        verbose = verbose,
    )
    ffs .*= factor
    return ffs #.* factor
end


function functionspace(basisfunctions)
    return basisfunctions.functionspace
end
