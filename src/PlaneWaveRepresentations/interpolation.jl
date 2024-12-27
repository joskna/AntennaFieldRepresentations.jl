abstract type InterpolateMap{Y<:SphereSamplingStrategy,T<:Real} <: LinearMaps.LinearMap{T} end

struct LocalθLocalΦInterpolateMap{Y<:SphereSamplingStrategy,orderθ,orderϕ,T<:Real} <:
       InterpolateMap{Y,T}
    originalsamplingstrategy::Y
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{SVector{orderθ,Int}}
    ϕweights::Vector{SVector{orderϕ,T}}
    ϕindices::Vector{SVector{orderϕ,Int}}
    ϕindicesopposite::Vector{SVector{orderϕ,Int}}
    posθranges::Vector{UnitRange{Int}}
    negθranges::Vector{UnitRange{Int}}
    intermediatestorage::Vector{Complex{T}}
    directions::Vector{Tuple{T,T}}
end
function Base.size(im::LocalθLocalΦInterpolateMap)
    θs, ϕs = samples(im.originalsamplingstrategy)
    return 2 * length(im.directions), 2 * length(θs) * length(ϕs)
end
function LocalθLocalΦInterpolateMap(
    originalsamplingstrategy::Y,
    directions::Vector{Tuple{T,T}};
    orderθ = 12,
    orderϕ = 12,
) where {Y<:SphereSamplingStrategy,T<:Real}
    return LocalθLocalΦInterpolateMap{Float64}(
        originalsamplingstrategy;
        directions,
        orderθ = orderθ,
        orderϕ = orderϕ,
    )
end
function LocalθLocalΦInterpolateMap{T}(
    originalsamplingstrategy::Y,
    directions::Vector{Tuple{T,T}};
    orderθ = 12,
    orderϕ = 12,
) where {T<:Real,Y<:SphereSamplingStrategy}
    oldθs, oldϕs = samples(originalsamplingstrategy)
    newθs = [direction[1] for direction in directions]
    newϕs = [direction[2] for direction in directions]


    intermediatestorage = zeros(Complex{T}, orderθ)
    # finalstorage = zeros(Complex{T}, length(directions))

    θweights, θindices, posθranges, negθranges =
        _planθweightsandindices(newθs, oldθs, orderθ, T)

    ϕweights, ϕindices = _planϕweightsandindices(newϕs, oldϕs, orderϕ, T)
    ϕindicesopposite::Vector{Svector{orderϕ,Int}}(undef, length(ϕindices))
    Nϕ = length(oldϕs)
    for (k, iϕrange) in enumerate(ϕindices)
        ϕindicesopposite[k] = SVector{orderϕ}(map_ϕrange!(iϕrange .+ (Nϕ ÷ 2), Nϕ))
    end

    return LocalθLocalΦInterpolateMap{Y,orderθ,orderϕ,T}(
        originalsamplingstrategy,
        targetsamplingstrategy,
        θweights,
        θindices,
        ϕweights,
        ϕindices,
        ϕindicesopposite,
        posθranges,
        negθranges,
        intermediatestorage,
        directions,
    )
end

abstract type ResampleMap{Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,T<:Real} <:
              LinearMaps.LinearMap{T} end

struct LocalθLocalΦResampleMap{
    Y1<:SphereSamplingStrategy,
    Y2<:SphereSamplingStrategy,
    orderθ,
    orderϕ,
    T<:Real,
} <: ResampleMap{Y1,Y2,T}
    originalsamplingstrategy::Y1
    targetsamplingstrategy::Y2
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{SVector{orderθ,Int}}
    ϕweights::Vector{SVector{orderϕ,T}}
    ϕindices::Vector{SVector{orderϕ,Int}}
    posθranges::Vector{UnitRange{Int}}
    negθranges::Vector{UnitRange{Int}}
    intermediatestorage::Matrix{Complex{T}}
    finalstorage::Matrix{Complex{T}}
end
function Base.size(rsm::LocalθLocalΦResampleMap)
    θsold, ϕsold = samples(rsm.originalsamplingstrategy)
    θsnew, ϕsnew = samples(rsm.targetsamplingstrategy)
    return 2 * length(θsnew) * length(ϕsnew), 2 * length(θsold) * length(ϕsold)
end
function LocalθLocalΦResampleMap(
    targetsamplingstrategy::Y2,
    originalsamplingstrategy::Y1;
    orderθ = 12,
    orderϕ = 12,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy}
    return LocalθLocalΦResampleMap{Float64}(
        targetsamplingstrategy,
        originalsamplingstrategy;
        orderθ = orderθ,
        orderϕ = orderϕ,
    )
end
function LocalθLocalΦResampleMap{T}(
    targetsamplingstrategy::Y2,
    originalsamplingstrategy::Y1;
    orderθ = 12,
    orderϕ = 12,
) where {T<:Real,Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy}
    oldθs, oldϕs = samples(originalsamplingstrategy)
    newθs, newϕs = samples(targetsamplingstrategy)


    intermediatestorage = zeros(Complex{T}, length(oldθs), length(newϕs))
    finalstorage = zeros(Complex{T}, length(newθs), length(newϕs))

    θweights, θindices, posθranges, negθranges =
        _planθweightsandindices(newθs, oldθs, orderθ, T)

    ϕweights, ϕindices = _planϕweightsandindices(newϕs, oldϕs, orderϕ, T)

    return LocalθLocalΦResampleMap{Y1,Y2,orderθ,orderϕ,T}(
        originalsamplingstrategy,
        targetsamplingstrategy,
        θweights,
        θindices,
        ϕweights,
        ϕindices,
        posθranges,
        negθranges,
        intermediatestorage,
        finalstorage,
    )
end

function _planϕweightsandindices(newϕs, oldϕs, orderϕ, T)
    ϕweights = Vector{SVector{orderϕ,T}}(undef, length(newϕs))
    ϕindices = Vector{SVector{orderϕ,Int64}}(undef, length(newϕs))
    Δϕ = oldϕs[2] - oldϕs[1]
    bϕ = lagrange_barycentric_weights((collect(1:orderϕ) .- 1) * Δϕ)
    for (k, newϕ) in enumerate(newϕs)
        iϕ₀ = find_next_smaller_ϕind(Δϕ, newϕ)
        iϕrange = ((iϕ₀-orderϕ+1):iϕ₀) .+ div(orderϕ, 2)

        ϕweights[k] = SVector{orderϕ}(
            T.(
                lagrange_interpolation_weights_from_barycentric(
                    (collect(iϕrange) .- 1) * Δϕ,
                    newϕ,
                    bϕ,
                )
            ),
        )
        ϕindices[k] = SVector{orderϕ}(map_ϕrange!(iϕrange, length(oldϕs)))

    end
    return ϕweights, ϕindices
end

function _planθweightsandindices(θvec_new, θvec_old, orderθ, T)
    newθlength = length(θvec_new)
    nθ_old = length(θvec_old)
    θweights = Vector{SVector{orderθ,T}}(undef, newθlength)
    θindices = Vector{MVector{orderθ,Int64}}(undef, newθlength)
    wθ_pos = Vector{SVector}(undef, newθlength)
    wθ_neg = Vector{SVector}(undef, newθlength)
    θinds_pos = Vector{SVector}(undef, newθlength)
    θinds_neg = Vector{SVector}(undef, newθlength)
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

    θweights = Vector{SVector{orderθ,T}}(undef, length(θvec_new))
    θindices = Vector{SVector{orderθ,Int64}}(undef, length(θvec_new))
    posθranges = Vector{UnitRange{Int}}(undef, length(θvec_new))
    negθranges = Vector{UnitRange{Int}}(undef, length(θvec_new))
    for k in eachindex(θindices)
        posθranges[k] = 1:length(θinds_pos[k])
        negθranges[k] = length(θinds_pos[k])+1:orderθ
        θindices[k] = SVector{orderθ,Int64}([θinds_pos[k]; θinds_neg[k]])
        θweights[k] = SVector{orderθ,T}([wθ_pos[k]; wθ_neg[k]])
    end

    return θweights, θindices, posθranges, negθranges
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

function extract_single_entry!(
    storage,
    Ematr,
    iθrange,
    posθrange,
    negθrange,
    wθ,
    iϕrange,
    ϕindicesopposite,
    wϕ,
)

    θinds_pos = view(iθrange, posθrange)
    θinds_neg = view(iθrange, negθrange)

    wθ_pos = view(wθ, posθrange)
    wθ_neg = view(wθ, negθrange)

    @inbounds mul!(storage, transpose(wθ_pos), view(Ematr, θinds_pos, iϕrange))
    @inbounds mul!(
        storage,
        transpose(wθ_neg),
        view(Ematr, θinds_neg, ϕindicesopposite),
        1,
        1,
    )
    return udot(wϕ, storage)
end

"""
    extract_single_planewave(Eθmatr::AbstractMatrix{<:Complex}, Eϕmatr::AbstractMatrix{<:Complex}, iθrange, wθ, iϕrange, wϕ ) -> Eθ, Eϕ

Use interpolation weights `wθ`, `wϕ` and correctly mapped ranges 'iθrange', `iϕrange` to extract Eθ- and Eϕ- value for single plane wave from `PlaneWaveExpansion`.
"""
function extract_single_planewave(
    Eθmatr::AbstractMatrix{C},
    Eϕmatr::AbstractMatrix{C},
    iθrange,
    wθ::AbstractArray{T},
    iϕrange,
    wϕ::AbstractArray{T},
) where {C<:Complex,T<:Real}

    Nθ, Nϕ = size(Eθmatr)
    iθrange = map_θrange!(iθrange, Nθ)
    iϕrange = map_ϕrange!(iϕrange, Nϕ)

    Eθvals = Vector{C}(undef, length(iθrange))
    Eϕvals = Vector{C}(undef, length(iθrange))
    iϕrange_new = map_ϕrange!(iϕrange .+ (Nϕ ÷ 2), Nϕ)

    for (k, θind) in enumerate(iθrange)
        if θind > 0
            Eθvals[k] = udot(wϕ, view(Eθmatr, θind, iϕrange))
            Eϕvals[k] = udot(wϕ, view(Eϕmatr, θind, iϕrange))
        else
            Eθvals[k] = -udot(wϕ, view(Eθmatr, -θind + 1, iϕrange_new))
            Eϕvals[k] = -udot(wϕ, view(Eϕmatr, -θind + 1, iϕrange_new))
        end
    end
    return udot(wθ, Eθvals), udot(wθ, Eϕvals)
end

# """
#     _extract_θ_row!(storage, Ematr::AbstractMatrix{C}, iθrange, wθ::<:AbstractVector) where {C<:Complex}

# Use interpolation weights `wθ` and correctly mapped range 'iθrange' to extract values for all plane waves in certain θ-cut from matrix of `PlaneWaveExpansion`.
# """
# function _extract_θ_row!(
#     Eθvals,
#     Eϕvals,
#     Eθmatr::AbstractMatrix{C},
#     Eϕmatr::AbstractMatrix{C},
#     iθrange,
#     wθ::AbstractVector,
# ) where {C<:Complex}

#     Nθ, Nϕ = size(Eθmatr)
#     iθrange = map_θrange!(iθrange, Nθ)
#     iϕrange = 1:Nϕ
#     iϕrange_new = map_ϕrange!((iϕrange) .+ (Nϕ ÷ 2), Nϕ)

#     ks = [k for (k, θind) in enumerate(iθrange) if θind > 0]
#     θinds = [θind for θind in iθrange if θind > 0]
#     # Eθvals .= transpose(wθ[ks]) * view(Eθmatr,θinds, :)
#     # Eϕvals .= transpose(wθ[ks]) * view(Eϕmatr,θinds, :)
#     mul!(Eθvals, transpose(wθ[ks]), view(Eθmatr, θinds, :))
#     mul!(Eϕvals, transpose(wθ[ks]), view(Eϕmatr, θinds, :))

#     ks = [k for (k, θind) in enumerate(iθrange) if θind <= 0]
#     θinds = [θind for θind in iθrange if θind <= 0]
#     Eθvals .+= -transpose(wθ[ks]) * view(Eθmatr, -θinds .+ 1, iϕrange_new)
#     Eϕvals .+= -transpose(wθ[ks]) * view(Eϕmatr, -θinds .+ 1, iϕrange_new)

#     return Eθvals, Eϕvals
# end


function _extract_single_θ!(
    storage::AbstractMatrix,
    Ematr::AbstractMatrix,
    iθrange,
    posθrange,
    negθrange,
    wθ,
)# where {C<:Complex}

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

function _extract_single_ϕ!(storage, Ematr, iϕrange, wϕ)
    mul!(storage, view(Ematr, :, iϕrange), wϕ)

    return storage
end
function _adjoint_extract_single_ϕ!(storage, Ematr, iϕrange, wϕ)
    view(Ematr, :, iϕrange) .+= storage .* adjoint(wϕ)
end

function _local_interpolate_θ!(
    Enew,
    E_old,
    newθlength,
    θindices,
    posθranges,
    negθranges,
    θweights,
)

    for k = 1:newθlength
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
function _adjoint_local_interpolate_θ!(
    Enew,
    E_old,
    newθlength,
    θindices,
    posθranges,
    negθranges,
    θweights;
    reset::Bool = true,
)
    if reset
        E_old .= zero(eltype(E_old))
    end
    for k = 1:newθlength
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


function _local_interpolate_ϕ!(Enew, E_old, newϕlength, ϕindices, ϕweights)
    for k = 1:newϕlength
        _extract_single_ϕ!(view(Enew, :, k), E_old, ϕindices[k], ϕweights[k])
    end
    return Enew
end
function _adjoint_local_interpolate_ϕ!(Enew, E_old, newϕlength, ϕindices, ϕweights)
    E_old .= zero(eltype(E_old))
    for k = 1:newϕlength
        _adjoint_extract_single_ϕ!(view(Enew, :, k), E_old, ϕindices[k], ϕweights[k])
    end
    return E_old
end

function _resamplematrix!(oldmatrix, interpolator::ResampleMap)
    return _resamplematrix!(interpolator.finalstorage, oldmatrix, interpolator)
end
function _resamplematrix!(
    storage,
    oldmatrix,
    interpolator::LocalθLocalΦResampleMap{Y1,Y2,orderθ,orderϕ,T},
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,orderθ,orderϕ,T}
    newθlength, newϕlength = size(interpolator.finalstorage)
    _local_interpolate_θ!(
        interpolator.intermediatestorage,
        oldmatrix,
        newθlength,
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
    )
    _local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage,
        newϕlength,
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    return storage
end
function _adjoint_resamplematrix!(
    storage,
    oldmatrix,
    interpolator::LocalθLocalΦResampleMap{Y1,Y2,orderθ,orderϕ,T};
    reset::Bool = true,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,orderθ,orderϕ,T}
    newθlength, newϕlength = size(interpolator.finalstorage)
    interpolator.intermediatestorage .= _adjoint_local_interpolate_ϕ!(
        storage,
        interpolator.intermediatestorage,
        newθlength,
        interpolator.ϕindices,
        interpolator.ϕweights,
    )
    oldmatrix .= _adjoint_local_interpolate_θ!(
        interpolator.intermediatestorage,
        oldmatrix,
        newϕlength,
        interpolator.θindices,
        interpolator.posθranges,
        interpolator.negθranges,
        interpolator.θweights,
        reset = reset,
    )
end

#Assume that interpolation coefficients are real valued
function _transpose_resamplematrix!(storage, oldmatrix, interpolator; reset = true)
    _adjoint_resamplematrix!(storage, oldmatrix, interpolator; reset = reset)
end

"""
    interpolate_single_planewave( (θnew,ϕnew), pattern::PlaneWaveExpansion, ::AbstractSphereInterpolation)

Interpolate single plane wave into direction `θnew`, `ϕnew` from given `PlaneWaveExpansion` using plynomial (Legendre) interpolation of order `orderθ` in θ and `orderϕ` in ϕ. 
"""
function interpolate_single_planewave(
    θnewϕnew::Tuple{T,T},
    pattern::PlaneWaveExpansion,
    ::Type{LocalθLocalΦInterpolateMap{Y,orderθ,orderϕ,T}};
    θvecϕvec = samples(pattern.samplingstrategy),
) where {Y,orderθ,orderϕ,T}
    θnew, ϕnew = θnewϕnew
    θvec, ϕvec = θvecϕvec

    Δϕ = ϕvec[2] - ϕvec[1]

    iθ₀ = find_next_smaller_θind(θvec, θnew) + orderθ ÷ 2
    iϕ₀ = find_next_smaller_ϕind(Δϕ, ϕnew) + orderϕ ÷ 2

    iθrange = ((iθ₀-orderθ+1):iθ₀)
    iϕrange = ((iϕ₀-orderϕ+1):iϕ₀)

    wϕ = lagrange_interpolation_weights(((iϕrange) .- 1) * Δϕ, ϕnew)
    wθ = lagrange_interpolation_weights(θsequence(θvec, iθrange), θnew)

    return extract_single_planewave(_eθ(pattern), _eϕ(pattern), iθrange, wθ, iϕrange, wϕ)

end

function interpolate_planewaves(
    directions::AbstractVector{Tuple{T,T}},
    pattern::PlaneWaveExpansion,
    ::Type{LocalθLocalΦInterpolateMap};
    orderθ = 12,
    orderϕ = 12,
) where {T}
    im = InterpolateMap(θϕs, pattern, orderθ = orderθ, orderϕ = orderϕ)
    FθFϕs = Matrix{eltype{pattern.EθEϕ}}(undef, length(directions), 2)
    FθFϕs .= interpolate_planewaves!(FθFϕs, im, pattern.EθEϕ)
    return [(FθFϕs[k, 1], FθFϕs[k, 2]) for k in eachindex(directions)]
end
function interpolate_planewaves!(
    FθFϕs,
    im::LocalθLocalΦInterpolateMap,
    EθEϕ::AbstractArray{C,3},
) where {C}
    for k in eachindex(im.directions)
        FθFϕs[k, 1] = extract_single_entry!(
            im.intermediatestorage,
            view(EθEϕ, :, :, 1),
            im.θindices[k],
            im.posθranges[k],
            im.negθranges[k],
            θweights[k],
            im.ϕindices[k],
            im.ϕindicesopposite[k],
            ϕweights[k],
        )
        FθFϕs[k, 2] = extract_single_entry!(
            im.intermediatestorage,
            view(EθEϕ, :, :, 2),
            im.θindices[k],
            im.posθranges[k],
            im.negθranges[k],
            θweights[k],
            im.ϕindices[k],
            im.ϕindicesopposite[k],
            ϕweights[k],
        )
    end
    return FθFϕs
end

function LinearMaps._unsafe_mul!(y, im::LocalθLocalΦInterpolateMap, x::AbstractVector)
    θs, ϕs = samples(im.originalsamplingstrategy)
    mat = reshape(x, length(θs), length(ϕs), 2)
    FθFϕs = reshape(y, length(im.directions), 2)
    FθFϕs .= interpolate_planewaves!(FθFϕs, im, mat)
    return y
end
