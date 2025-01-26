
abstract type ResampleMap{Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,T<:Real} <:
              LinearMaps.LinearMap{T} end

struct LocalθResampleMap{
    Y1<:SphereSamplingStrategy,
    Y2<:SphereSamplingStrategy,
    orderθ,
    T<:Real,
} <: ResampleMap{Y1,Y2,T}
    originalsamplingstrategy::Y1
    targetsamplingstrategy::Y2
    θweights::Vector{SVector{orderθ,T}}
    θindices::Vector{SVector{orderθ,Int}}
    posθranges::Vector{UnitRange{Int}}
    negθranges::Vector{UnitRange{Int}}
    inputbuffer::Vector{Complex{T}}
    outputbuffer::Vector{Complex{T}}
    inputbuffermat::Array{Complex{T},3}
    outputbuffermat::Array{Complex{T},3}
    storage::SubArray{Complex{T},2}
end
function LocalθResampleMap(
    targetsamplingstrategy::Y2,
    originalsamplingstrategy::Y1;
    orderθ = 12,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy}
    T = Float64
    oldθs, oldϕs = samples(originalsamplingstrategy)
    newθs, newϕs = samples(targetsamplingstrategy)

    @assert norm(oldϕs - newϕs) < 1e-3

    inputbuffer = zeros(Complex{T}, 2 * length(oldθs) * length(oldϕs))
    outputbuffer = zeros(Complex{T}, 2 * length(newθs) * length(newϕs))
    inputbuffermat = reshape(inputbuffer, length(oldθs), length(oldϕs), 2)
    outputbuffermat = reshape(outputbuffer, length(newθs), length(newϕs), 2)
    storage = view(outputbuffermat, :, :, 1)
    θweights, θindices, posθranges, negθranges =
        _planθweightsandindices(newθs, oldθs, orderθ, T)

    return LocalθResampleMap{Y1,Y2,orderθ,T}(
        originalsamplingstrategy,
        targetsamplingstrategy,
        θweights,
        θindices,
        posθranges,
        negθranges,
        inputbuffer,
        outputbuffer,
        inputbuffermat,
        outputbuffermat,
        storage,
    )
end

struct LocalϕResampleMap{
    Y1<:SphereSamplingStrategy,
    Y2<:SphereSamplingStrategy,
    orderϕ,
    T<:Real,
} <: ResampleMap{Y1,Y2,T}
    originalsamplingstrategy::Y1
    targetsamplingstrategy::Y2
    ϕweights::Vector{SVector{orderϕ,T}}
    ϕindices::Vector{SVector{orderϕ,Int}}
    inputbuffer::Vector{Complex{T}}
    outputbuffer::Vector{Complex{T}}
    inputbuffermat::Array{Complex{T},3}
    outputbuffermat::Array{Complex{T},3}
    storage::SubArray{Complex{T},2}
end
function LocalϕResampleMap(
    targetsamplingstrategy::Y2,
    originalsamplingstrategy::Y1;
    orderϕ = 12,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy}
    T = Float64
    oldθs, oldϕs = samples(originalsamplingstrategy)
    newθs, newϕs = samples(targetsamplingstrategy)

    @assert norm(oldθs - newθs) < 1e-3

    inputbuffer = zeros(Complex{T}, 2 * length(oldθs) * length(oldϕs))
    outputbuffer = zeros(Complex{T}, 2 * length(newθs) * length(newϕs))
    inputbuffermat = reshape(inputbuffer, length(oldθs), length(oldϕs), 2)
    outputbuffermat = reshape(outputbuffer, length(newθs), length(newϕs), 2)
    storage = view(outputbuffermat, :, :, 1)
    ϕweights, ϕindices = _planϕweightsandindices(newϕs, oldϕs, orderϕ, T)

    return LocalϕResampleMap{Y1,Y2,orderϕ,T}(
        originalsamplingstrategy,
        targetsamplingstrategy,
        ϕweights,
        ϕindices,
        inputbuffer,
        outputbuffer,
        inputbuffermat,
        outputbuffermat,
        storage,
    )
end

struct θϕResampleMap{
    R1<:ResampleMap,
    R2<:ResampleMap,
    Y1<:SphereSamplingStrategy,
    Y2<:SphereSamplingStrategy,
    orderθ,
    orderϕ,
    T<:Real,
} <: ResampleMap{Y1,Y2,T}
    originalsamplingstrategy::Y1
    targetsamplingstrategy::Y2
    θresamplemap::R1
    ϕresamplemap::R2
    inputbuffer::Vector{Complex{T}}
    outputbuffer::Vector{Complex{T}}
    inputbuffermat::Array{Complex{T},3}
    outputbuffermat::Array{Complex{T},3}
end

LocalθLocalϕResampleMap{Y1,Y2,orderθ,orderϕ,T} =
    θϕResampleMap{LocalθResampleMap,LocalϕResampleMap,Y1,Y2,orderθ,orderϕ,T}

function Base.size(rsm::ResampleMap)
    return length(rsm.outputbuffer), length(rsm.inputbuffer)
end

function LocalθLocalϕResampleMap(
    targetsamplingstrategy::Y2,
    originalsamplingstrategy::Y1;
    orderθ = 12,
    orderϕ = 12,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy}

    T = Float64
    intermediatesamplingstrategy =
        _intermediate_samplingstrategy(targetsamplingstrategy, originalsamplingstrategy)

    θresamplemap = LocalθResampleMap(
        intermediatesamplingstrategy,
        originalsamplingstrategy;
        orderθ = orderθ,
    )

    ϕresamplemap = LocalϕResampleMap(
        targetsamplingstrategy,
        intermediatesamplingstrategy;
        orderϕ = orderϕ,
    )

    inputbuffer = θresamplemap.inputbuffer
    outputbuffer = ϕresamplemap.outputbuffer
    inputbuffermat = θresamplemap.inputbuffermat
    outputbuffermat = ϕresamplemap.outputbuffermat



    return LocalθLocalϕResampleMap{Y1,Y2,orderθ,orderϕ,T}(
        originalsamplingstrategy,
        targetsamplingstrategy,
        θresamplemap,
        ϕresamplemap,
        inputbuffer,
        outputbuffer,
        inputbuffermat,
        outputbuffermat,
    )
end


function _intermediate_samplingstrategy(
    targetsamplingstrategy::GaussLegendreθRegularϕSampling,
    originalsamplingstrategy::GaussLegendreθRegularϕSampling,
)
    Jϕ_old = originalsamplingstrategy.Jϕ

    Nθ_new = targetsamplingstrategy.Nθ

    return GaussLegendreθRegularϕSampling(Nθ_new, Jϕ_old)
end

function _intermediate_samplingstrategy(
    targetsamplingstrategy::RegularθRegularϕSampling,
    originalsamplingstrategy::GaussLegendreθRegularϕSampling,
)
    Jϕ_old = originalsamplingstrategy.Jϕ

    Jθ_new = targetsamplingstrategy.Jθ

    return RegularθRegularϕSampling(Jθ_new, Jϕ_old)
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
    rsm::LocalθLocalϕResampleMap{Y1,Y2,orderθ,orderϕ,T},
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,orderθ,orderϕ,T}
    newθlength, newϕlength = size(interpolator.finalstorage)
    _local_interpolate_θ!(
        rsm.θresamplemap.storage,
        oldmatrix,
        newθlength,
        rsm.θindices,
        rsm.posθranges,
        rsm.negθranges,
        rsm.θweights,
    )
    _local_interpolate_ϕ!(
        storage,
        rsm.θresamplemap.storage,
        newϕlength,
        rsm.ϕindices,
        rsm.ϕweights,
    )
    return storage
end
function _adjoint_resamplematrix!(
    storage,
    oldmatrix,
    rsm::LocalθLocalϕResampleMap{Y1,Y2,orderθ,orderϕ,T};
    reset::Bool = true,
) where {Y1<:SphereSamplingStrategy,Y2<:SphereSamplingStrategy,orderθ,orderϕ,T}
    newθlength, newϕlength = size(rsm.finalstorage)
    rsm.θresamplemap.storage .= _adjoint_local_interpolate_ϕ!(
        storage,
        rsm.θresamplemap.storage,
        newθlength,
        rsm.ϕindices,
        rsm.ϕweights,
    )
    oldmatrix .= _adjoint_local_interpolate_θ!(
        rsm.θresamplemap.storage,
        oldmatrix,
        newϕlength,
        rsm.θindices,
        rsm.posθranges,
        rsm.negθranges,
        rsm.θweights,
        reset = reset,
    )
end

#Assume that interpolation coefficients are real valued
function _transpose_resamplematrix!(storage, oldmatrix, rsm; reset = true)
    _adjoint_resamplematrix!(storage, oldmatrix, rsm; reset = reset)
end

function LinearMaps._unsafe_mul!(y, rsm::LocalθResampleMap, x::AbstractVector)

    # θs, ϕs = samples(rsm.originalsamplingstrategy)
    # newθs, newϕs = samples(rsm.targetsamplingstrategy)
    # mat = reshape(x, length(θs), length(ϕs), 2)


    # matout = reshape(y, length(newθs), length(newϕs), 2)

    # setindex! to rsm.inputbuffer also writes into rsm.inputbuffermat
    rsm.inputbuffer .= x

    # setindex! to rsm.outputbuffer also writes into rsm.outputbuffermat
    rsm.outputbuffer .= 0

    _local_interpolate_θ!(
        view(rsm.outputbuffermat, :, :, 1),
        view(rsm.inputbuffermat, :, :, 1),
        size(rsm.outputbuffermat, 1),
        rsm.θindices,
        rsm.posθranges,
        rsm.negθranges,
        rsm.θweights,
    )

    _local_interpolate_θ!(
        view(rsm.outputbuffermat, :, :, 2),
        view(rsm.inputbuffermat, :, :, 2),
        size(rsm.outputbuffermat, 1),
        rsm.θindices,
        rsm.posθranges,
        rsm.negθranges,
        rsm.θweights,
    )

    y .= rsm.outputbuffer
    return y
end


function LinearMaps._unsafe_mul!(y, rsm::LocalϕResampleMap, x::AbstractVector)
    # θs, ϕs = samples(rsm.originalsamplingstrategy)
    # newθs, newϕs = samples(rsm.targetsamplingstrategy)
    # mat = reshape(x, length(θs), length(ϕs), 2)

    # # setindex! to matout writes into y
    # matout = reshape(y, length(newθs), length(newϕs), 2)

    # rsm.storage .= 0

    # setindex! to rsm.inputbuffer also writes into rsm.inputbuffermat
    rsm.inputbuffer .= x

    # setindex! to rsm.outputbuffer also writes into rsm.outputbuffermat
    rsm.outputbuffer .= 0
    _local_interpolate_ϕ!(
        view(rsm.outputbuffermat, :, :, 1),
        view(rsm.inputbuffermat, :, :, 1),
        size(rsm.outputbuffermat, 2),
        rsm.ϕindices,
        rsm.ϕweights,
    )

    _local_interpolate_ϕ!(
        view(rsm.outputbuffermat, :, :, 2),
        view(rsm.inputbuffermat, :, :, 2),
        size(rsm.outputbuffermat, 2),
        rsm.ϕindices,
        rsm.ϕweights,
    )
    y .= rsm.outputbuffer


    return y
end

function LinearMaps._unsafe_mul!(y, rsm::θϕResampleMap, x::AbstractVector)

    rsm.inputbuffer .= x

    rsm.θresamplemap.outputbuffer .= rsm.θresamplemap * rsm.inputbuffer

    rsm.outputbuffer .= rsm.ϕresamplemap * rsm.θresamplemap.outputbuffer

    y .= rsm.outputbuffer

    return y
end
