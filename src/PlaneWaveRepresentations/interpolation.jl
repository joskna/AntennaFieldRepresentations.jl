abstract type AbstractInterpolationStrategy end

struct LocalInterpolation{Lold,Lnew,orderθ,orderϕ,T} <:
       AbstractInterpolationStrategy where {
    Lold<:Integer,
    Lnew<:Integer,
    orderθ<:Integer,
    orderϕ<:Integer,
    T<:Real,
} end
function Lold(
    ::LocalInterpolation{Lold_,Lnew,orderθ,orderϕ,T},
) where {Lold_,Lnew,orderθ,orderϕ,T}
    return Lold_
end
function Lnew(
    ::LocalInterpolation{Lold,Lnew_,orderθ,orderϕ,T},
) where {Lold,Lnew_,orderθ,orderϕ,T}
    return Lnew_
end

struct LocalThetaGlobalPhiInterpolation{Lold,Lnew,orderθ,T} <:
       AbstractInterpolationStrategy where {
    Lold<:Integer,
    Lnew<:Integer,
    orderθ<:Integer,
    T<:Real,
} end
function Lold(
    ::LocalThetaGlobalPhiInterpolation{Lold_,Lnew,orderθ,T},
) where {Lold_,Lnew,orderθ,T}
    return Lold_
end
function Lnew(
    ::LocalThetaGlobalPhiInterpolation{Lold,Lnew_,orderθ,T},
) where {Lold,Lnew_,orderθ,T}
    return Lnew_
end

struct GlobalInterpolation{Lold,Lnew,T} <:
       AbstractInterpolationStrategy where {Lold<:Integer,Lnew<:Integer,T<:Real}
    finalstorage::Matrix{Complex{T}}
    sphericalexpansion::RadiatingSphericalExpansion{Complex{T}}
end
function GlobalInterpolation{T}(Lold, Lnew) where {T<:Real}
    finalstorage = Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2)
    sphericalexpansion = RadiatingSphericalExpansion{Complex{T}}(
        Vector{Complex{T}}(undef, sℓm_to_j(2, Lold, Lold)),
    )
    return GlobalInterpolation{Lold,Lnew,T}(finalstorage, sphericalexpansion)
end
function Lold(::GlobalInterpolation{Lold_,Lnew,T}) where {Lold_,Lnew,T}
    return Lold_
end
function Lnew(::GlobalInterpolation{Lold,Lnew_,T}) where {Lold,Lnew_,T}
    return Lnew_
end
struct PlannedLocalInterpolation{orderθ,orderϕ,T} <:
       AbstractInterpolationStrategy where {orderθ<:Integer,orderϕ<:Integer,T<:Real}
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{SVector{orderθ,Int}}
    ϕweights::Vector{SVector{orderϕ,T}}
    ϕindices::Vector{SVector{orderϕ,Int}}
    posθranges::Vector{UnitRange{Int}}
    negθranges::Vector{UnitRange{Int}}
    intermediatestorage::Matrix{Complex{T}}
    finalstorage::Matrix{Complex{T}}
    Lold::Int
    Lnew::Int
end
struct ParallelPlannedLocalInterpolation{orderθ,orderϕ,T} <:
       AbstractInterpolationStrategy where {orderθ<:Integer,orderϕ<:Integer,T<:Real}
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{SVector{orderθ,Int}}
    ϕweights::Vector{SVector{orderϕ,T}}
    ϕindices::Vector{SVector{orderϕ,Int}}
    posθranges::Vector{UnitRange{Int}}
    negθranges::Vector{UnitRange{Int}}
    intermediatestorage::Vector{Matrix{Complex{T}}}
    finalstorage::Vector{Matrix{Complex{T}}}
    Lold::Int
    Lnew::Int
end


function PlannedLocalInterpolation{T}(Lold, Lnew; orderθ = 12, orderϕ = 12) where {T<:Real}
    intermediatestorage = Array{Complex{T}}(undef, Lnew + 1, 2 * Lold + 2)
    finalstorage = Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2)
    # intermediatestorage=zeros(Complex{T}, Lnew+1, 2*Lold+2)
    θweights = Vector{SVector{orderθ,T}}(undef, Lnew + 1)
    θindices = Vector{SVector{orderθ,Int64}}(undef, Lnew + 1)
    ϕweights = Vector{SVector{orderϕ,T}}(undef, 2 * Lnew + 2)
    ϕindices = Vector{SVector{orderθ,Int64}}(undef, 2 * Lnew + 2)
    _, θvec_new, ϕvec_new = samplingrule(Lnew)
    _, θvec_old, ϕvec_old = samplingrule(Lold)
    # nθ_old = Lold + 1

    posθranges = Vector{UnitRange{Int}}(undef, Lnew + 1)
    negθranges = Vector{UnitRange{Int}}(undef, Lnew + 1)
    wθ_pos, wθ_neg, θinds_pos, θinds_neg =
        _planθweightsandindices(θvec_new, θvec_old, orderθ, T)
    for k in eachindex(θindices)
        posθranges[k] = 1:length(θinds_pos[k])
        # negθranges[k]=length(θinds_pos[k]) < orderθ ? (length(θinds_pos[k])+1 : orderθ) : (UnitRange{Int}[])
        negθranges[k] = (length(θinds_pos[k])+1:orderθ)
        θindices[k] = SVector{orderθ,Int64}([θinds_pos[k]; θinds_neg[k]])
        θweights[k] = SVector{orderθ,T}([wθ_pos[k]; wθ_neg[k]])
    end

    # indold = 0
    # maxindex = Int(ceil(length(θvec_new) / 2))
    # for kθ in 1:maxindex
    #     startindex = maximum([indold, 1])
    #     θnew = θvec_new[kθ]
    #     indold = find_next_smaller_θind(θvec_old[startindex:end], θnew) + startindex
    #     iθ₀ = indold + orderθ ÷ 2
    #     iθrange = map_θrange!(((iθ₀ - orderθ + 1):iθ₀), nθ_old)
    #     iθrange_mirrored = -(iθrange) .+ (nθ_old + 1)
    #     wθ = lagrange_interpolation_weights(θsequence(θvec_old, iθrange), θnew)
    #     # wθ= _sign_θinterpolationweights!(wθ, iθrange, nθ_old)

    #     θweights[kθ] = SVector{orderθ}(T.(wθ))
    #     θindices[kθ] = iθrange

    #     θweights[end - kθ + 1] = SVector{orderθ}(T.(wθ))
    #     θindices[end - kθ + 1] = iθrange_mirrored

    # end

    # _, __, ϕvec_old = samplingrule(Lold)
    Δϕ = ϕvec_old[2] - ϕvec_old[1]
    bϕ = lagrange_barycentric_weights((collect(1:orderϕ) .- 1) * Δϕ)
    for (k, ϕnew) in enumerate(ϕvec_new)
        iϕ₀ = find_next_smaller_ϕind(Δϕ, ϕnew)
        iϕrange = ((iϕ₀-orderϕ+1):iϕ₀) .+ div(orderϕ, 2)

        ϕweights[k] = SVector{orderϕ}(
            T.(
                lagrange_interpolation_weights_from_barycentric(
                    (collect(iϕrange) .- 1) * Δϕ,
                    ϕnew,
                    bϕ,
                )
            ),
        )
        ϕindices[k] = SVector{orderϕ}(map_ϕrange!(iϕrange, length(ϕvec_old)))

    end
    return PlannedLocalInterpolation{orderθ,orderϕ,T}(
        θweights,
        θindices,
        ϕweights,
        ϕindices,
        posθranges,
        negθranges,
        intermediatestorage,
        finalstorage,
        Lold,
        Lnew,
    )

end
function ParallelPlannedLocalInterpolation{T}(
    Lold,
    Lnew;
    orderθ = 12,
    orderϕ = 12,
) where {T<:Real}
    intermediatestorage =
        [Array{Complex{T}}(undef, Lnew + 1, 2 * Lold + 2) for k = 1:Threads.nthreads()]
    finalstorage =
        [Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2) for k = 1:Threads.nthreads()]
    # intermediatestorage=zeros(Complex{T}, Lnew+1, 2*Lold+2)
    θweights = Vector{SVector{orderθ,T}}(undef, Lnew + 1)
    θindices = Vector{SVector{orderθ,Int64}}(undef, Lnew + 1)
    ϕweights = Vector{SVector{orderϕ,T}}(undef, 2 * Lnew + 2)
    ϕindices = Vector{SVector{orderθ,Int64}}(undef, 2 * Lnew + 2)
    _, θvec_new, ϕvec_new = samplingrule(Lnew)
    _, θvec_old, ϕvec_old = samplingrule(Lold)
    # nθ_old = Lold + 1

    posθranges = Vector{UnitRange{Int}}(undef, Lnew + 1)
    negθranges = Vector{UnitRange{Int}}(undef, Lnew + 1)
    wθ_pos, wθ_neg, θinds_pos, θinds_neg =
        _planθweightsandindices(θvec_new, θvec_old, orderθ, T)
    for k in eachindex(θindices)
        posθranges[k] = 1:length(θinds_pos[k])
        # negθranges[k]=length(θinds_pos[k]) < orderθ ? (length(θinds_pos[k])+1 : orderθ) : (UnitRange{Int}[])
        negθranges[k] = (length(θinds_pos[k])+1:orderθ)
        θindices[k] = SVector{orderθ,Int64}([θinds_pos[k]; θinds_neg[k]])
        θweights[k] = SVector{orderθ,T}([wθ_pos[k]; wθ_neg[k]])
    end
    Δϕ = ϕvec_old[2] - ϕvec_old[1]
    bϕ = lagrange_barycentric_weights((collect(1:orderϕ) .- 1) * Δϕ)
    for (k, ϕnew) in enumerate(ϕvec_new)
        iϕ₀ = find_next_smaller_ϕind(Δϕ, ϕnew)
        iϕrange = ((iϕ₀-orderϕ+1):iϕ₀) .+ div(orderϕ, 2)

        ϕweights[k] = SVector{orderϕ}(
            T.(
                lagrange_interpolation_weights_from_barycentric(
                    (collect(iϕrange) .- 1) * Δϕ,
                    ϕnew,
                    bϕ,
                )
            ),
        )
        ϕindices[k] = SVector{orderϕ}(map_ϕrange!(iϕrange, length(ϕvec_old)))

    end
    return ParallelPlannedLocalInterpolation{orderθ,orderϕ,T}(
        θweights,
        θindices,
        ϕweights,
        ϕindices,
        posθranges,
        negθranges,
        intermediatestorage,
        finalstorage,
        Lold,
        Lnew,
    )

end

import Base.show
function Base.show(
    io::IO,
    interpolator::PlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    print(io, "PlannedLocalInterpolation{$orderθ,$orderϕ,$T}")
end
function Lold(P::PlannedLocalInterpolation{orderθ,orderϕ,T}) where {orderθ,orderϕ,T}
    return P.Lold
end
function Lnew(P::PlannedLocalInterpolation{orderθ,orderϕ,T}) where {orderθ,orderϕ,T}
    return P.Lnew
end
function Base.show(
    io::IO,
    interpolator::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    print(io, "ParallelPlannedLocalInterpolation{$orderθ,$orderϕ,$T}")
end
function Lold(P::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T}) where {orderθ,orderϕ,T}
    return P.Lold
end
function Lnew(P::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T}) where {orderθ,orderϕ,T}
    return P.Lnew
end

struct PlannedLocalThetaGlobalPhiInterpolation{Lold,Lnew,orderθ,T} <:
       AbstractInterpolationStrategy where {
    Lold<:Integer,
    Lnew<:Integer,
    orderθ<:Integer,
    T<:Real,
}
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{MVector{orderθ,Int64}}
    fftplan!::AbstractFFTs.Plan
    ifftplan!::AbstractFFTs.Plan
    intermediatestorage::Matrix{Complex{T}}
    finalstorage::Matrix{Complex{T}}
end
function PlannedLocalThetaGlobalPhiInterpolation{T}(Lold, Lnew; orderθ = 12) where {T<:Real}
    intermediatestorage = Array{Complex{T}}(undef, Lnew + 1, 2 * Lold + 2)
    finalstorage = Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2)
    # θweights = Vector{SVector{orderθ,T}}(undef, Lnew + 1)
    # θindices = Vector{MVector{orderθ,Int64}}(undef, Lnew + 1)
    _, θvec_new, __ = samplingrule(Lnew)
    _, θvec_old, __ = samplingrule(Lold)
    # nθ_old = Lold + 1

    # indold = 0
    # maxindex = Int(ceil(length(θvec_new) / 2))
    # for kθ in 1:maxindex
    #     startindex = maximum([indold, 1])
    #     θnew = θvec_new[kθ]
    #     indold = find_next_smaller_θind(θvec_old[startindex:end], θnew) + startindex
    #     iθ₀ = indold + orderθ ÷ 2
    #     iθrange = map_θrange!(((iθ₀ - orderθ + 1):iθ₀), nθ_old)
    #     iθrange_mirrored = -(iθrange) .+ (nθ_old + 1)
    #     wθ = lagrange_interpolation_weights(θsequence(θvec_old, iθrange), θnew)
    #     # wθ= _sign_θinterpolationweights!(wθ, iθrange, nθ_old)

    #     θweights[kθ] = SVector{orderθ,T}(wθ)
    #     θindices[kθ] = iθrange

    #     θweights[end - kθ + 1] = SVector{orderθ}(convert.(T, wθ))
    #     θindices[end - kθ + 1] = iθrange_mirrored
    # end
    θweights, θindices = _planθweightsandindices(θvec_new, θvec_old, orderθ, T)

    fftplan! = plan_fft!(intermediatestorage, 2)
    ifftplan! = plan_ifft!(Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2), 2)
    return PlannedLocalThetaGlobalPhiInterpolation{Lold,Lnew,orderθ,T}(
        θweights,
        θindices,
        fftplan!,
        ifftplan!,
        intermediatestorage,
        finalstorage,
    )
end
function _planθweightsandindices(θvec_new, θvec_old, orderθ, T)
    Lnew = length(θvec_new) - 1
    nθ_old = length(θvec_old)
    θweights = Vector{SVector{orderθ,T}}(undef, Lnew + 1)
    θindices = Vector{MVector{orderθ,Int64}}(undef, Lnew + 1)
    wθ_pos = Vector{SVector}(undef, Lnew + 1)
    wθ_neg = Vector{SVector}(undef, Lnew + 1)
    θinds_pos = Vector{SVector}(undef, Lnew + 1)
    θinds_neg = Vector{SVector}(undef, Lnew + 1)
    indold = 0
    maxindex = Int(ceil(length(θvec_new) / 2))
    for kθ = 1:maxindex
        startindex = maximum([indold, 1])
        θnew = θvec_new[kθ]
        indold = find_next_smaller_θind(θvec_old[startindex:end], θnew) + startindex
        iθ₀ = indold + orderθ ÷ 2
        # iθrange = map_θrange!(((iθ₀ - orderθ + 1):iθ₀), nθ_old)
        # iθrange_mirrored = -(iθrange) .+ (nθ_old + 1)
        iθrange = ((iθ₀-orderθ+1):iθ₀)
        wθ = lagrange_interpolation_weights(θsequence(θvec_old, iθrange), θnew)
        iθrange = map_θrange!(((iθ₀-orderθ+1):iθ₀), nθ_old)
        iθrange_mirrored = -(iθrange) .+ (nθ_old + 1)
        # wθ= _sign_θinterpolationweights!(wθ, iθrange, nθ_old)

        θweights[kθ] = SVector{orderθ,T}(wθ)
        θindices[kθ] = iθrange

        θweights[end-kθ+1] = SVector{orderθ}(convert.(T, wθ))
        θindices[end-kθ+1] = iθrange_mirrored

    end
    for kθ in eachindex(θvec_new)
        wθ_pos_tmp, wθ_neg_tmp, θinds_pos_tmp, θinds_neg_tmp =
            _split_θ_interpolationparams(θweights[kθ], θindices[kθ], nθ_old)
        wθ_pos[kθ] = SVector{length(wθ_pos_tmp)}(wθ_pos_tmp)
        wθ_neg[kθ] = SVector{length(wθ_neg_tmp)}(wθ_neg_tmp)
        θinds_pos[kθ] = SVector{length(wθ_pos_tmp)}(θinds_pos_tmp)
        θinds_neg[kθ] = SVector{length(wθ_neg_tmp)}(θinds_neg_tmp)
    end
    return wθ_pos, wθ_neg, θinds_pos, θinds_neg
end
function Lold(
    ::PlannedLocalThetaGlobalPhiInterpolation{Lold_,Lnew,orderθ,T},
) where {Lold_,Lnew,orderθ,T}
    return Lold_
end
function Lnew(
    ::PlannedLocalThetaGlobalPhiInterpolation{Lold,Lnew_,orderθ,T},
) where {Lold,Lnew_,orderθ,T}
    return Lnew_
end

function _sign_θinterpolationweights!(wθ, iθrange, Nθ)
    mappediθrange = deepcopy(iθrange)
    mappediθrange = map_θrange!(mappediθrange, Nθ)
    ks = [k for (k, θind) in enumerate(mappediθrange) if θind <= 0]
    wθ[ks] .*= -1
    return wθ
end



"""
    lagrange_interpolation_weights(X::Array{::Number,1},Y::Array{::Number,1},t::Number) -> w::Array{Flota64}

Return interpolation weights wᵢ for the calculation f(t)= ∑ᴺᵢ₌₁ wᵢ⋅ f(xᵢ) 

# Extended help 
The interpolation weights are obtained using the "second (true) form of the barycentric formula".
See J.-P. Berrut and L. N. Trefethen: "Barycentric Lagrange Interpolation", SIAM Review, Vol. 46, No. 3, pp. 501-517
"""
function lagrange_interpolation_weights(X, t::T) where {T<:Number}
    b = lagrange_barycentric_weights(X)
    return lagrange_interpolation_weights_from_barycentric(X, t, b)
end
function lagrange_interpolation_weights_from_barycentric(
    X,
    t::T,
    barycentric_weights::Vector{<:Number},
) where {T<:Number}
    w = zeros(T, length(barycentric_weights))
    for (k, b) in enumerate(barycentric_weights)
        if t ≈ X[k]
            w = zeros(T, length(barycentric_weights))
            w[k] = 1
            return w
        end
        w[k] = b / (t - X[k])
    end
    return w / sum(w)
end

"""
    lagrange_barycentric_weights(X::Array{<:Number,1}) -> w::Array{<:Number,1}

Return barycentric weights for calculating Lagrange interpolation weights

# Extended help 
The barycentric weights are needed to evaluate the "second (true) form of the barycentric formula" for the Lagrange interpolation.
See J.-P. Berrut and L. N. Trefethen: "Barycentric Lagrange Interpolation", SIAM Review, Vol. 46, No. 3, pp. 501-517
"""
function lagrange_barycentric_weights(X)
    T = eltype(X)
    b = zeros(T, length(X))
    b[1] = 1
    for j = 2:length(X)
        for k = 1:(j-1)
            b[k] = (X[k] - X[j]) * b[k]
        end
        b[j] = prod((X[j] - X[k]) for k = 1:(j-1))
    end
    return 1 ./ b
end

"""
    find_next_smaller_θind(θvec::Array{<:Real,1}, θnew<:Real) -> i::Integer

Return index of largest entry in ascendingly ordered θvec which is still smaller than θnew
"""
function find_next_smaller_θind(θvec::Array{<:Real,1}, θnew::Real)
    return searchsortedlast(θvec, θnew)
end

"""
    find_next_smaller_ϕind(Δϕ<:Real, ϕnew<:Real)-> i::Integer

Return index of largest entry in ascendingly ordered ϕvec which is still smaller than ϕnew
"""
function find_next_smaller_ϕind(Δϕ::Real, ϕnew::Real)
    return convert(Int64, floor(ϕnew / Δϕ))
end


"""
    θsequence(θvec::Array{<:Real, 1}, iθrange::Array{<:Integer,1})

Truncate and/or periodically extend `θvec` to create `θsequence` matching the index range in `iθrange`.
"""
function θsequence(θvec::Array{<:Real,1}, iθrange)
    if iθrange[1] < 1
        # attach shifted copy of θvec before θvec
        # maximum shiftvalue must be smaller than min(θvec)
        shiftvalue = ((θvec[end] - θvec[1]) ÷ pi + 1) * pi
        newθvec = [θvec .- shiftvalue; θvec]
        newiθrange = iθrange .+ length(θvec)
        return θsequence(newθvec, newiθrange)
    elseif iθrange[end] > length(θvec)
        # attach shifted copy of θvec before θvec
        # maximum shiftvalue must be smaller than min(θvec)
        shiftvalue = ((θvec[end] - θvec[1]) ÷ pi + 1) * pi
        newθvec = [θvec; θvec .+ shiftvalue]
        return θsequence(newθvec, iθrange)
    end
    # Nθ=length(θvec)
    # mappediθrange=map_θrange!(copy(iθrange), Nθ)
    # ks = [k for (k, θind) in enumerate(mappediθrange) if θind <= 0]
    # mappediθrange[ks]= -mappediθrange[ks] .+ 1

    return θvec[iθrange]
end


"""
    map_θrange!(iθrange, Nθ)

Map the index range `iθrange` to values between `-Nθ+1` and  `Nθ`.
"""
function map_θrange!(iθrange, Nθ::Integer)
    iθrange = mod.(iθrange, 2 * Nθ)
    for k in eachindex(iθrange)
        if iθrange[k] > Nθ
            iθrange[k] += -2 * Nθ
        elseif iθrange[k] < -Nθ + 1
            iθrange[k] += 2 * Nθ
        end
    end
    return iθrange
end

"""
    map_ϕrange!(iϕrange, Nθ)

Map the index range `iϕrange` to values between `1` and  `Nϕ`.
"""
function map_ϕrange!(iϕrange, Nϕ::Integer)
    iϕrange_vector = mod.(iϕrange, Nϕ)
    for (kϕ, iϕ) in enumerate(iϕrange_vector)
        if iϕ == 0
            iϕrange_vector[kϕ] = Nϕ
        end
    end
    return iϕrange_vector
end


"""
    extract_single_planewave(Eθmatr::Array{<:Complex,2}, Eϕmatr::Array{<:Complex,2}, iθrange, wθ, iϕrange, wϕ ) -> Eθ, Eϕ

Use interpolation weights `wθ`, `wϕ` and correctly mapped ranges 'iθrange', `iϕrange` to extract Eθ- and Eϕ- value for single plane wave from `PlaneWaveRepresentation`.
"""
function extract_single_planewave(
    Eθmatr::Array{C,2},
    Eϕmatr::Array{C,2},
    iθrange,
    wθ::Array{<:Real,1},
    iϕrange,
    wϕ::Array{<:Real,1},
) where {C<:Complex}

    Nθ, Nϕ = size(Eθmatr)
    iθrange = map_θrange!(iθrange, Nθ)
    iϕrange = map_ϕrange!(iϕrange, Nϕ)

    Eθvals = Vector{C}(undef, length(iθrange))
    Eϕvals = Vector{C}(undef, length(iθrange))
    iϕrange_new = map_ϕrange!(iϕrange .+ (Nϕ ÷ 2), Nϕ)

    for (k, θind) in enumerate(iθrange)
        if θind > 0
            Eθvals[k] = udot(wϕ, Eθmatr[θind, iϕrange])
            Eϕvals[k] = udot(wϕ, Eϕmatr[θind, iϕrange])
        else
            Eθvals[k] = -udot(wϕ, Eθmatr[-θind+1, iϕrange_new])
            Eϕvals[k] = -udot(wϕ, Eϕmatr[-θind+1, iϕrange_new])
        end
    end
    return udot(wθ, Eθvals), udot(wθ, Eϕvals)
end

"""
    _extract_θ_row!(storage, Ematr::Array{C,2}, iθrange, wθ::<:AbstractVector) where {C<:Complex}

Use interpolation weights `wθ` and correctly mapped range 'iθrange' to extract values for all plane waves in certain θ-cut from matrix of `PlaneWaveRepresentation`.
"""
function _extract_θ_row!(
    Eθvals,
    Eϕvals,
    Eθmatr::Array{C,2},
    Eϕmatr::Array{C,2},
    iθrange,
    wθ::AbstractVector,
) where {C<:Complex}

    Nθ, Nϕ = size(Eθmatr)
    iθrange = map_θrange!(iθrange, Nθ)
    iϕrange = 1:Nϕ
    iϕrange_new = map_ϕrange!((iϕrange) .+ (Nϕ ÷ 2), Nϕ)

    ks = [k for (k, θind) in enumerate(iθrange) if θind > 0]
    θinds = [θind for θind in iθrange if θind > 0]
    # Eθvals .= transpose(wθ[ks]) * view(Eθmatr,θinds, :)
    # Eϕvals .= transpose(wθ[ks]) * view(Eϕmatr,θinds, :)
    mul!(Eθvals, transpose(wθ[ks]), view(Eθmatr, θinds, :))
    mul!(Eϕvals, transpose(wθ[ks]), view(Eϕmatr, θinds, :))

    ks = [k for (k, θind) in enumerate(iθrange) if θind <= 0]
    θinds = [θind for θind in iθrange if θind <= 0]
    Eθvals .+= -transpose(wθ[ks]) * view(Eθmatr, -θinds .+ 1, iϕrange_new)
    Eϕvals .+= -transpose(wθ[ks]) * view(Eϕmatr, -θinds .+ 1, iϕrange_new)

    return Eθvals, Eϕvals
end

function _split_θ_interpolationparams(wθ, iθrange, Nθ)
    # Nθ, Nϕ = size(Ematr)
    iθrange .= map_θrange!(iθrange, Nθ)

    ks = [k for (k, θind) in enumerate(iθrange) if θind > 0]
    ks2 = [k for (k, θind) in enumerate(iθrange) if θind <= 0]

    θinds_pos = [θind for θind in iθrange if θind > 0]
    θinds_neg = -[θind for θind in iθrange if θind <= 0] .+ 1

    wθ_pos = view(wθ, ks)
    wθ_neg = -view(wθ, ks2)

    return wθ_pos, wθ_neg, θinds_pos, θinds_neg
end

function _iϕranges(Nϕ)
    iϕrange_pos = 1:Nϕ
    Nϕhalf = Nϕ ÷ 2
    iϕrange_neg = [Nϕhalf+1:Nϕ; 1:Nϕhalf]
    # iϕrange_neg = map_ϕrange!((iϕrange_pos) .+ (Nϕ ÷ 2), Nϕ)
    return iϕrange_pos, iϕrange_neg
end


function _extract_single_θ!(storage, Ematr, iθrange, posθrange, negθrange, wθ)# where {C<:Complex}

    Nϕ = size(Ematr, 2)
    Nϕhalf = Nϕ ÷ 2

    θinds_pos = view(iθrange, posθrange)
    θinds_neg = view(iθrange, negθrange)

    wθ_pos = view(wθ, posθrange)
    wθ_neg = view(wθ, negθrange)

    @inbounds mul!(storage, transpose(wθ_pos), view(Ematr, θinds_pos, 1:Nϕ))
    @inbounds mul!(
        view(storage, :, 1:Nϕhalf),
        transpose(wθ_neg),
        view(Ematr, θinds_neg, Nϕhalf+1:Nϕ),
        1,
        1,
    )
    @inbounds mul!(
        view(storage, :, Nϕhalf+1:Nϕ),
        transpose(wθ_neg),
        view(Ematr, θinds_neg, 1:Nϕhalf),
        1,
        1,
    )

    return storage
end
# function _extract_single_θ!(storage, Ematr::Array{C,2}, θinds_pos, θinds_neg, wθ_pos, wθ_neg, iϕrange_pos, iϕrange_neg) where {C<:Complex}

#     mul!(storage, transpose(wθ_pos), view(Ematr, θinds_pos, iϕrange_pos))
#     mul!(storage, transpose(wθ_neg), view(Ematr, θinds_neg, iϕrange_neg), 1,1)    

#     return storage
# end


function _adjoint_extract_single_θ!(storage, Ematr, iθrange, posθrange, negθrange, wθ)
    Nϕ = size(Ematr, 2)
    Nϕhalf = Nϕ ÷ 2

    θinds_pos = view(iθrange, posθrange)
    θinds_neg = view(iθrange, negθrange)

    wθ_pos = view(wθ, posθrange)
    wθ_neg = view(wθ, negθrange)

    @inbounds for k in eachindex(θinds_pos), kk = 1:Nϕ
        Ematr[θinds_pos[k], kk] += wθ_pos[k] * storage[kk]
    end
    @inbounds for k in eachindex(θinds_neg), kk = 1:Nϕhalf
        Ematr[θinds_neg[k], kk] += wθ_neg[k] * storage[kk+Nϕhalf]
        Ematr[θinds_neg[k], kk+Nϕhalf] += wθ_neg[k] * storage[kk]
    end


end
# function _adjoint_extract_single_θ!(storage, Ematr::AbstractMatrix, θinds_pos, θinds_neg, wθ_pos, wθ_neg, iϕrange_pos, iϕrange_neg)

#     # view(Ematr,θinds_pos, iϕrange_pos) .+= conj.(wθ_pos) * transpose(storage)
#     # view(Ematr,θinds_neg, iϕrange_neg) .+= -conj.(wθ_neg) * transpose(storage)
#     # mul!(view(Ematr,θinds_pos, iϕrange_pos), conj.(wθ_pos), transpose(storage), 1 ,1)
#     # mul!(view(Ematr,θinds_neg, iϕrange_neg), conj.(wθ_neg), transpose(storage), -1 ,1)

#     view(Ematr,θinds_pos, iϕrange_pos) .+= conj.(wθ_pos) * transpose(storage)
#     view(Ematr,θinds_neg, iϕrange_neg) .+= conj.(wθ_neg) * transpose(storage)
#     # view(Ematr,θinds_neg, iϕrange_neg) .+=  transpose(storage * transpose(conj.(wθ_neg)))
#     # mul!(view(Ematr,θinds_neg, iϕrange_neg), conj.(wθ_neg), transpose(storage), 1 ,1)
#     # mul!(transpose(view(Ematr,θinds_neg, iϕrange_neg)), storage, transpose(conj.(wθ_neg)), 1 ,1)

# end
function _extract_single_ϕ!(storage, Ematr, iϕrange, wϕ)
    # _, Nϕ = size(Ematr)
    # iϕrange .= map_ϕrange!(iϕrange, Nϕ)

    # storage .= Ematr[:, map_ϕrange!(copy(iϕrange), Nϕ)] * wϕ
    # mul!(storage, view(Ematr, :, map_ϕrange!(iϕrange, Nϕ)), wϕ)
    # storage .= view(Ematr, :, iϕrange) * wϕ
    mul!(storage, view(Ematr, :, iϕrange), wϕ)

    return storage
end
function _adjoint_extract_single_ϕ!(storage, Ematr, iϕrange, wϕ)
    view(Ematr, :, iϕrange) .+= storage .* adjoint(wϕ)
    # mul!(view(Ematr, :, iϕrange) , storage , adjoint(wϕ))
end


function _local_interpolate_ϕ(
    Eθ_old,
    Eϕ_old,
    ϕvec_old::AbstractVector,
    ϕvec_new::AbstractVector,
    interpolationorder::Integer,
)
    @assert size(Eθ_old) == size(Eϕ_old)
    nθ_old, nϕ_old = size(Eθ_old)
    @assert nϕ_old == length(ϕvec_old)

    Eθ_new = similar(Eθ_old, nθ_old, length(ϕvec_new))
    Eϕ_new = similar(Eϕ_old, nθ_old, length(ϕvec_new))

    Δϕ = ϕvec_old[2] - ϕvec_old[1]
    bϕ = lagrange_barycentric_weights((collect(1:interpolationorder) .- 1) * Δϕ)
    for (k, ϕnew) in enumerate(ϕvec_new)
        iϕ₀ = find_next_smaller_ϕind(Δϕ, ϕnew) + interpolationorder ÷ 2
        iϕrange = ((iϕ₀-interpolationorder+1):iϕ₀)

        wϕ = lagrange_interpolation_weights_from_barycentric(
            (collect(iϕrange) .- 1) * Δϕ,
            ϕnew,
            bϕ,
        )
        iϕrange = map_ϕrange!(iϕrange, length(ϕvec_old))

        Eθ_new[:, k] .= view(Eθ_old, :, iϕrange) * wϕ
        Eϕ_new[:, k] .= view(Eϕ_old, :, iϕrange) * wϕ
    end
    return Eθ_new, Eϕ_new
end

function _global_interpolate_ϕ(E_old, nϕ_old::Integer, nϕ_new::Integer)

    @assert nϕ_old == size(E_old, 2)

    # Eθ_new = resample(Eθ_old, (nθ_old, length(ϕvec_new)))
    # Eϕ_new = resample(Eϕ_old, (nθ_old, length(ϕvec_new)))
    E_new = _resampleϕ_fft(E_old, nϕ_new)

    return E_new
end



function _resampleϕ_fft(Aold, nΦnew)
    nθold, nΦold = size(Aold)

    if nΦnew > nΦold
        _upsampleϕ_fft!(Aold, nΦnew, nθold, nΦold)
    elseif nΦnew < nΦold
        _downsampleϕ_fft!(Aold, nΦnew, nθold, nΦold)
    else
        return Aold
    end
end
function _upsampleϕ_fft!(Aold, nΦnew, nθold, nΦold)
    Afft = fft!(Aold, 2)
    Afftnew = fill!(similar(Aold, nθold, nΦnew), 0)
    # Afft=view(Afftnew,:,[1:div(nΦold,2); nΦold-cld(nΦold,2)+1:nΦold])
    # Afft.=fft!(Aold,2)        
    Afftnew[:, 1:div(nΦold, 2)] .= Afft[:, 1:div(nΦold, 2)]
    Afftnew[:, (nΦnew-cld(nΦold, 2)+1):nΦnew] .= Afft[:, (nΦold-cld(nΦold, 2)+1):nΦold]
    return nΦnew / nΦold .* ifft!(Afftnew, 2)
end
function _downsampleϕ_fft!(Aold, nΦnew, nθold, nΦold)
    Afft = fft!(Aold, 2)
    Afftnew = fill!(similar(Aold, nθold, nΦnew), 0)
    Afftnew[:, 1:div(nΦnew, 2)] .= Afft[:, 1:div(nΦnew, 2)]
    Afftnew[:, (nΦnew-cld(nΦnew, 2)+1):nΦnew] .= Afft[:, (nΦold-cld(nΦnew, 2)+1):nΦold]

    # Afftnew=view(Afft, [1:div(nΦnew,2); nΦold-cld(nΦnew,2)+1:nΦold])
    return nΦnew / nΦold .* ifft!(Afftnew, 2)
end

function _planned_resample_fft!(Afftnew, Aold, nΦnew, fftplan!, ifftplan!)
    _, nΦold = size(Aold)

    if nΦnew > nΦold
        _planned_upsampleϕ_fft!(Afftnew, Aold, nΦnew, nΦold, fftplan!, ifftplan!)
    elseif nΦnew < nΦold
        _planned_downsampleϕ_fft!(Afftnew, Aold, nΦnew, nΦold, fftplan!, ifftplan!)
    else
        return Aold
    end
end
function _planned_upsampleϕ_fft!(Afftnew, Aold, nΦnew, nΦold, fftplan!, ifftplan!)
    Afft = fftplan! * Aold
    fill!(Afftnew, 0)
    # Afft=view(Afftnew,:,[1:div(nΦold,2); nΦold-cld(nΦold,2)+1:nΦold])
    # Afft.=fft!(Aold,2)        
    Afftnew[:, 1:div(nΦold, 2)] .= view(Afft, :, 1:div(nΦold, 2))
    Afftnew[:, (nΦnew-cld(nΦold, 2)+1):nΦnew] .=
        view(Afft, :, (nΦold-cld(nΦold, 2)+1):nΦold)
    Afftnew .= ifftplan! * Afftnew
    Afftnew .*= (nΦnew / nΦold)
    return Afftnew
end
function _planned_downsampleϕ_fft!(Afftnew, Aold, nΦnew, nΦold, fftplan!, ifftplan!)
    Afft = fftplan! * Aold
    fill!(Afftnew, 0)
    # Afft=view(Afftnew,:,[1:div(nΦold,2); nΦold-cld(nΦold,2)+1:nΦold])
    # Afft.=fft!(Aold,2)        
    Afftnew[:, 1:div(nΦnew, 2)] .= view(Afft, :, 1:div(nΦnew, 2))
    Afftnew[:, (nΦnew-cld(nΦnew, 2)+1):nΦnew] .=
        view(Afft, :, (nΦold-cld(nΦnew, 2)+1):nΦold)
    Afftnew .= ifftplan! * Afftnew
    Afftnew .*= (nΦnew / nΦold)
    return Afftnew
end


function _local_interpolate_θ(
    Eθ_old,
    Eϕ_old,
    θvec_old::AbstractVector,
    θvec_new::AbstractVector,
    interpolationorder::Integer,
)
    nθ_old, nϕ_old = size(Eθ_old)
    Eθ = similar(Eθ_old, length(θvec_new), nϕ_old)
    Eϕ = similar(Eϕ_old, length(θvec_new), nϕ_old)
    indold = 0
    maxindex = Int(ceil(length(θvec_new) / 2))
    Eθvals, Eϕvals = similar(Eθ_old, 1, nϕ_old), similar(Eϕ_old, 1, nϕ_old)
    for kθ = 1:maxindex
        startindex = maximum([indold, 1])
        θnew = θvec_new[kθ]
        indold = find_next_smaller_θind(θvec_old[startindex:end], θnew) + startindex
        iθ₀ = indold + interpolationorder ÷ 2
        iθrange = map_θrange!(((iθ₀-interpolationorder+1):iθ₀), nθ_old)
        iθrange_mirrored = -(iθrange) .+ (nθ_old + 1)
        wθ = lagrange_interpolation_weights(θsequence(θvec_old, iθrange), θnew)

        Eθ[kθ, :], Eϕ[kθ, :] =
            _extract_θ_row!(Eθvals, Eϕvals, Eθ_old, Eϕ_old, (iθrange), wθ)
        Eθ[end-kθ+1, :], Eϕ[end-kθ+1, :] =
            _extract_θ_row!(Eθvals, Eϕvals, Eθ_old, Eϕ_old, iθrange_mirrored, wθ)
    end
    return Eθ, Eϕ
end

function _local_interpolate_θ!(
    Enew,
    E_old,
    Lnew,
    θindices,
    posθranges,
    negθranges,
    θweights,
)
    nθ_new = Lnew + 1

    for k = 1:nθ_new
        _extract_single_θ!(
            transpose(view(Enew, k, :)),
            E_old,
            θindices[k],
            posθranges[k],
            negθranges[k],
            θweights[k],
        )
    end
    return Enew
end

# function _local_interpolate_θ!(Enew, E_old, θinds_pos, θinds_neg, wθ_pos, wθ_neg)     
#     nθ_new, Nϕ = size(Enew)
#     iϕrange_pos, iϕrange_neg = _iϕranges(Nϕ)

#         for k in 1:nθ_new
#             _extract_single_θ!(transpose(view(Enew, k, :)), E_old, θinds_pos[k], θinds_neg[k], wθ_pos[k], wθ_neg[k], iϕrange_pos, iϕrange_neg)
#         end
#         return Enew
#     end

function _adjoint_local_interpolate_θ!(
    Enew,
    E_old,
    Lnew,
    θindices,
    posθranges,
    negθranges,
    θweights;
    reset::Bool = true,
)
    nθ_new = Lnew + 1
    if reset
        E_old .= zero(eltype(E_old))
    end
    for k = 1:nθ_new
        _adjoint_extract_single_θ!(
            view(Enew, k, :),
            E_old,
            θindices[k],
            posθranges[k],
            negθranges[k],
            θweights[k],
        )
    end
    return E_old
end
# function _adjoint_local_interpolate_θ!(Enew, E_old, Lnew, θinds_pos, θinds_neg, wθ_pos, wθ_neg; reset::Bool=true) 
#     nθ_new, Nϕ = size(Enew)
#     iϕrange_pos, iϕrange_neg = _iϕranges(Nϕ)
#     if reset 
#         E_old .= zero(eltype(E_old))
#     end
#     for k in 1:nθ_new
#         _adjoint_extract_single_θ!(view(Enew, k, :), E_old, θinds_pos[k], θinds_neg[k], wθ_pos[k], wθ_neg[k], iϕrange_pos, iϕrange_neg)
#     end
#     return E_old
# end

function _local_interpolate_ϕ!(Enew, E_old, Lnew, ϕindices, ϕweights)
    nϕ_new = 2 * Lnew + 2

    # Threads.@threads for chunk in collect(Iterators.partition(1:nϕ_new, nϕ_new ÷ Threads.nthreads() ))
    #     for k in chunk
    #         _extract_single_ϕ!(view(Enew, :, k), E_old, ϕindices[k], ϕweights[k])
    #     end
    # end
    for k = 1:nϕ_new
        _extract_single_ϕ!(view(Enew, :, k), E_old, ϕindices[k], ϕweights[k])
    end
    return Enew
end
function _adjoint_local_interpolate_ϕ!(Enew, E_old, Lnew, ϕindices, ϕweights)
    nϕ_new = 2 * Lnew + 2
    E_old .= zero(eltype(E_old))
    for k = 1:nϕ_new
        _adjoint_extract_single_ϕ!(view(Enew, :, k), E_old, ϕindices[k], ϕweights[k])
    end
    return E_old
end

function _interpolatematrix!(
    storage,
    oldmatrix,
    interpolator::PlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    _local_interpolate_θ!(
        interpolator.intermediatestorage,
        oldmatrix,
        Lnew(interpolator),
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
    )
    _local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage,
        Lnew(interpolator),
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    return storage
end
function _interpolatematrix!(
    oldmatrix,
    interpolator::PlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    _interpolatematrix!(interpolator.finalstorage, oldmatrix, interpolator)

end
function _interpolatematrix!(
    storage,
    oldmatrix,
    interpolator::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    _local_interpolate_θ!(
        interpolator.intermediatestorage[Threads.threadid()],
        oldmatrix,
        Lnew(interpolator),
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
    )
    _local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage[Threads.threadid()],
        Lnew(interpolator),
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    return storage
end
function _interpolatematrix!(
    oldmatrix,
    interpolator::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T},
) where {orderθ,orderϕ,T}
    _interpolatematrix!(
        interpolator.finalstorage[Threads.threadid()],
        oldmatrix,
        interpolator,
    )
end
function _interpolatematrix!(
    storage,
    oldmatrix,
    interpolator::PlannedLocalThetaGlobalPhiInterpolation{orderθ,T},
) where {orderθ,T}
    _local_interpolate_θ!(
        interpolator.intermediatestorage,
        oldmatrix,
        Lnew(interpolator),
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
    )
    _planned_resample_fft!(
        storage,
        interpolator.intermediatestorage,
        2 * Lnew(interpolator) + 2,
        interpolator.fftplan!,
        interpolator.ifftplan!,
    )
    return storage
end
function _interpolatematrix!(
    storage,
    oldmatrix,
    ::GlobalInterpolation{Lold,Lnew,T},
) where {Lold,Lnew,T}
    pattern = FarfieldPattern{Complex{T}}(Lold, oldmatrix, oldmatrix)
    α = convertrepresentation(RadiatingSphericalExpansion{Complex{T}}, pattern)
    newpattern = _convertrepresentation_and_resample(FarfieldPattern{Complex{T}}, α, Lnew)
    storage .= (newpattern.Eθ .+ newpattern.Eϕ) ./ 2
    return storage
end

function _adjoint_interpolatematrix!(
    storage,
    oldmatrix,
    interpolator::PlannedLocalInterpolation{orderθ,orderϕ,T};
    reset::Bool = true,
) where {orderθ,orderϕ,T}
    interpolator.intermediatestorage .= _adjoint_local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage,
        Lnew(interpolator),
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    oldmatrix .= _adjoint_local_interpolate_θ!(
        interpolator.intermediatestorage,
        oldmatrix,
        Lnew(interpolator),
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
        reset = reset,
    )
end
function _adjoint_interpolatematrix!(
    storage,
    oldmatrix,
    interpolator::ParallelPlannedLocalInterpolation{orderθ,orderϕ,T};
    reset::Bool = true,
) where {orderθ,orderϕ,T}
    interpolator.intermediatestorage[Threads.threadid()] .= _adjoint_local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage[Threads.threadid()],
        Lnew(interpolator),
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    oldmatrix .= _adjoint_local_interpolate_θ!(
        interpolator.intermediatestorage[Threads.threadid()],
        oldmatrix,
        Lnew(interpolator),
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
        reset = reset,
    )
end

#Assume that interpolation coefficients are real valued
function _transpose_interpolatematrix!(storage, oldmatrix, interpolator; reset = true)
    _adjoint_interpolatematrix!(storage, oldmatrix, interpolator; reset = reset)
end

"""
    resample(pattern::PlaneWaveRepresentation, L::Integer, [orderθ::Integer=12, orderϕ::Integer=12]) -> newpattern::PlaneWaveRepresentation

Return pattern sampled according to order `L` using Lagrange interpolation of orders `orderθ`, `orderϕ`     
"""
function resample(pattern::FarfieldPattern, interpolator::AbstractInterpolationStrategy)
    # @assert L == Lnew
    @assert pattern.L == Lold(interpolator)
    Eθnew = Array{Complex{T}}(undef, Lnew(interpolator) + 1, 2 * Lnew(interpolator) + 2)
    Eϕnew = Array{Complex{T}}(undef, Lnew(interpolator) + 1, 2 * Lnew(interpolator) + 2)

    Eθnew = _interpolatematrix!(Eθnew, pattern.Eθ, interpolator)
    Eϕnew = _interpolatematrix!(Eϕnew, pattern.Eϕ, interpolator)

    newpattern = FarfieldPattern(L, Eθnew, Eϕnew)
    # newpattern.L = Lnew
    # newpattern.Eθ .= Eθnew
    # newpattern.Eϕ .= Eϕnew

    return newpattern

end
# function resample(
#     pattern::FarfieldPattern, L::Integer, interpolator::PlannedLocalThetaGlobalPhiInterpolation{Lold,Lnew,orderθ,T}
# ) where {Lold,Lnew,orderθ,T}
#     @assert L == Lnew
#     @assert pattern.L == Lold
#     Eθnew = Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2)
#     Eϕnew = Array{Complex{T}}(undef, Lnew + 1, 2 * Lnew + 2)

#     _local_interpolate_θ!(interpolator.intermediatestorage, pattern.Eθ, Lnew, interpolator.θindices, interpolator.θweights)
#     _planned_resample_fft!(Eθnew, interpolator.intermediatestorage, 2 * Lnew + 2, interpolator.fftplan!, interpolator.ifftplan!)

#     _local_interpolate_θ!(interpolator.intermediatestorage, pattern.Eϕ, Lnew, interpolator.θindices, interpolator.θweights)
#     _planned_resample_fft!(Eϕnew, interpolator.intermediatestorage, 2 * Lnew + 2, interpolator.fftplan!, interpolator.ifftplan!)

#     newpattern = FarfieldPattern(L, Eθnew, Eϕnew)
#     # newpattern.L = Lnew
#     # newpattern.Eθ .= Eθnew
#     # newpattern.Eϕ .= Eϕnew

#     return newpattern

# end
function resample(
    pattern::FarfieldPattern,
    ::LocalInterpolation{Lold,Lnew,orderθ,orderϕ,T},
) where {Lold,Lnew,orderθ,orderϕ,T}
    @assert pattern.L == Lold
    newpattern = deepcopy(pattern)
    newpattern.L = Lnew
    _, θvec_old, ϕvec_old = samplingrule(Lold)
    _, θvec_new, ϕvec_new = samplingrule(Lnew)

    Eθaux, Eϕaux = _local_interpolate_θ(pattern.Eθ, pattern.Eϕ, θvec_old, θvec_new, orderθ)
    Eθ, Eϕ = _local_interpolate_ϕ(Eθaux, Eϕaux, ϕvec_old, ϕvec_new, orderϕ)

    newpattern.Eθ = Eθ
    newpattern.Eϕ = Eϕ
    return newpattern
end
function resample(
    pattern::PlaneWaveSpectrum,
    ::LocalInterpolation{Lold,Lnew,orderθ,orderϕ,T},
) where {Lold,Lnew,orderθ,orderϕ,T}
    @assert pattern.L == Lold
    newpattern = deepcopy(pattern)
    newpattern.L = Lnew
    _, θvec_old, ϕvec_old = samplingrule(Lold)
    _, θvec_new, ϕvec_new = samplingrule(Lnew)

    Eθaux, Eϕaux = _local_interpolate_θ(pattern.Eθ, pattern.Eϕ, θvec_old, θvec_new, orderθ)
    Eθ, Eϕ = _local_interpolate_ϕ(Eθaux, Eϕaux, ϕvec_old, ϕvec_new, orderϕ)

    newpattern.Eθ = Eθ
    newpattern.Eϕ = Eϕ
    return newpattern
end
function resample(pattern::PlaneWaveRepresentation, Lnew::Integer)
    return resample(
        pattern,
        LocalInterpolation{pattern.L,Lnew,12,12,real(elementtype(pattern))},
    )
end
function resample(
    pattern::PlaneWaveRepresentation,
    ::LocalThetaGlobalPhiInterpolation{Lold,Lnew,orderθ,T},
) where {Lold,Lnew,orderθ,T}
    # newpattern = deepcopy(pattern)
    # newpattern.L = L
    _, θvec_old, ϕvec_old = samplingrule(Lold)
    _, θvec_new, ϕvec_new = samplingrule(Lnew)

    Eθaux, Eϕaux = _local_interpolate_θ(pattern.Eθ, pattern.Eϕ, θvec_old, θvec_new, orderθ)
    Eθ = _global_interpolate_ϕ(Eθaux, length(ϕvec_old), length(ϕvec_new))
    Eϕ = _global_interpolate_ϕ(Eϕaux, length(ϕvec_old), length(ϕvec_new))

    return typeof(pattern)(L, Eθ, Eϕ)
end

function resample(
    pattern::FarfieldPattern,
    interpolator::GlobalInterpolation{Lold,Lnew,T},
) where {Lold,Lnew,T}
    @assert pattern.L == Lold
    interpolator.sphericalexpansion .=
        convertrepresentation(RadiatingSphericalExpansion{Complex{T}}, pattern)
    return _convertrepresentation_and_resample(
        typeof(pattern),
        interpolator.sphericalexpansion,
        Lnew,
    )
end


"""
    interpolate_single_planewave(θnew::Real, ϕnew::Real, pattern::PlaneWaveRepresentation, ::AbstractInterpolationStrategy)

Interpolate single plane wave into direction `θnew`, `ϕnew` from given `PlaneWaveRepresentation` using plynomial (Legendre) interpolation of order `orderθ` in θ and `orderϕ` in ϕ. 
"""
function interpolate_single_planewave(
    θnew::Real,
    ϕnew::Real,
    pattern::PlaneWaveRepresentation,
    ::LocalInterpolation{Lold,Lnew,orderθ,orderϕ,T},
) where {Lold,Lnew,orderθ,orderϕ,T}

    @assert pattern.L == Lold

    _, θvec, ϕvec = samplingrule(Lold)

    Δϕ = ϕvec[2] - ϕvec[1]

    iθ₀ = find_next_smaller_θind(θvec, θnew) + orderθ ÷ 2
    iϕ₀ = find_next_smaller_ϕind(Δϕ, ϕnew) + orderϕ ÷ 2

    iθrange = ((iθ₀-orderθ+1):iθ₀)
    iϕrange = ((iϕ₀-orderϕ+1):iϕ₀)

    wϕ = lagrange_interpolation_weights(((iϕrange) .- 1) * Δϕ, ϕnew)
    wθ = lagrange_interpolation_weights(θsequence(θvec, iθrange), θnew)

    return extract_single_planewave(pattern.Eθ, pattern.Eϕ, iθrange, wθ, iϕrange, wϕ)

end
