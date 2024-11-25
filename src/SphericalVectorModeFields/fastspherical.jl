####################### auxilliary functions #################################################
#                         |
#                         v
####################################################################################################

function _Πfun(ℓ)
    if isodd(ℓ)
        return 0.0
    else
        return 2 / (1 - ℓ^2)
    end
end

function _convfun!(storage, P, B, Bext, L, fftplan!, ifftplan)
    # length(A) == 4 * L + 1 ? () : throw(DimensionMismatch("length(A) must be 4L+1"))
    lenB = length(B)
    fill!(Bext, zero(eltype(Bext)))
    # numzeros= lenB == 2 * L + 1 ? 2*L : 2*L -1

    # if lenB != 2 * L + 1  && lenB != 2*L+2
    #      throw(DimensionMismatch("length(B) must be 2L+1 or 2L+2"))
    # end
    view(Bext, 1:L+1) .= view(B, 1:L+1)
    lenBext = length(Bext)
    view(Bext, lenBext-(lenB-(L+2)):lenBext) .= view(B, L+2:lenB)
    mul!(Bext, fftplan!, Bext)
    Bext .= P .* Bext
    mul!(storage, ifftplan, Bext)
    return storage
    # return ifftplan*( P .* (fftplan * Bext))
    # return ifft(P .* fft(Bext))
end

# @memoize function _Δℓ_mμ(ℓ)
function _Δℓ_mμ(ℓ)
    # println(string("Storing value for ℓ= ", ℓ))
    return wignerd(ℓ, pi / 2)
end

# for ℓ in 1:100
#     @memoize _Δℓ_mμ(ℓ)
# end

# @memoize function _Δℓ_mμ(ℓ, ϑ)
function _Δℓ_mμ(ℓ, ϑ)
    # println(string("Storing value for ℓ= ", ℓ))
    return wignerd(ℓ, ϑ)
end

# function S12_to_v(S12,L)
function _χandφ_integral(
    S12::AbstractArray{C,3},
    L::Integer,
    Jθ::Integer,
    Jϕ::Integer,
) where {C<:Complex}
    v = zeros(C, 2 * L + 1, 2 * L + 1, 2)
    Lθ = (Jθ - 1) ÷ 2
    vview = view(v, 1:(Lθ+1), 1:Jϕ, :)
    ifft_planϕ = plan_ifft!(vview, 2)
    return _χandφ_integral!(v, S12, L, ifft_planϕ, Jθ, Jϕ)
end
function _φ_integral!(
    v::AbstractArray{C,3},
    vview,
    L::Integer,
    ifft_planϕ,
    Jθ::Integer,
    Jϕ::Integer,
) where {C<:Complex}

    Lθ = (Jθ - 1) ÷ 2
    # Lϕ = (Jϕ - 1) ÷ 2
    # vview = view(v, 1:(Jθ-Lθ), :, :)   
    # vview .= mul!(vview, ifft_planϕ, vview)
    mul!(vview, ifft_planϕ, vview)

    # add mirror images of theta modes to neg. theta indices
    lenv = size(v, 1)
    view(v, (lenv-Lθ+1):lenv, :, :) .= view(v, (Lθ+1):-1:2, :, :)

    # overwrite redundant ϕ-modes
    view(v, :, L+2:2L+1, :) .= view(v, :, Jϕ-L+1:Jϕ, :)

    μval = Int.(fftfreq(2 * L + 1, 2 * L + 1))
    for k = 1:(2*L+1)
        # signfac = (-1)^(μval[k] + 1)
        signfac = _negpow1(μval[k] + 1)
        # v[(Lθ+2):end, k, :] .*= signfac
        view(v, (Lθ+2):size(v, 1), k, :) .*= signfac
    end
    return v
end

function _χ_integral!(w::AbstractArray{C,3}, S12::AbstractArray{C,3}) where {C}
    # w[:, :, 1] .= view(S12, :, :, 1) + C(0, 1) .* view(S12, :, :, 2) # w_1(ϑ, ϕ)
    # w[:, :, 2] .= view(S12, :, :, 1) + C(0, -1) .* view(S12, :, :, 2) # w_-1 (ϑ, ϕ)

    mul!(view(w, :, :, 1), UniformScaling(C(0, 1)), view(S12, :, :, 2))
    view(w, :, :, 1) .+= view(S12, :, :, 1)# w_1(ϑ, ϕ)

    mul!(view(w, :, :, 2), UniformScaling(C(0, -1)), view(S12, :, :, 2))
    view(w, :, :, 2) .+= view(S12, :, :, 1)# w_-1(ϑ, ϕ)


    w .*= 0.5

    return w
end


# function v_to_u(v,L)
# function _θintegral(v::AbstractArray{C,3}, L::Integer) where {C<:Number}

#     u = zeros(C, L * (L + 2), 2)
#     P = _Πfun.(ifftshift((-2*L):(2*L)))
#     K = zeros(C, 4 * L + 1, 2 * L + 1, 2)

#     Δ = Vector{Matrix{Float64}}(undef, L)
#     for ℓ = 1:L
#         Δ[ℓ] = Matrix{Float64}(undef, 2ℓ + 1, 2ℓ + 1)
#         ifftshift!(Δ[ℓ], _Δℓ_mμ(ℓ))
#     end

#     return _θintegral!(u, v, L, P, K, Δ)
# end
function _θintegral!(
    u::AbstractArray{C},
    v::AbstractArray{C,3},
    L::Integer,
    P,
    K,
    Δ,
    ifft_planθ,
) where {C<:Number}
    #regularly distributed θ
    fill!(u, zero(C))
    Jθ = size(v, 1)

    mul!(v, ifft_planθ, v)
    vview = view(v, [1:L+1; Jθ-L+1:Jθ], :, :)

    Bext = zeros(C, size(P))
    fftplan! = plan_fft!(Bext)
    ifftplan = plan_ifft(Bext)
    for m_ind = 1:(2*L+1)
        _convfun!(
            view(K, :, m_ind, 1),
            P,
            view(vview, :, m_ind, 1),
            Bext,
            L,
            fftplan!,
            ifftplan,
        )
        _convfun!(
            view(K, :, m_ind, 2),
            P,
            view(vview, :, m_ind, 2),
            Bext,
            L,
            fftplan!,
            ifftplan,
        )

        # _convfun!(K[:, m_ind, 1], P, vview[:, m_ind, 1], Bext, L, fftplan, ifftplan)
        # _convfun!(K[:, m_ind, 2], P, vview[:, m_ind, 2], Bext, L, fftplan, ifftplan)
    end

    # K[2:end, :, :] .*= 2
    view(K, 2:size(K, 1), :, :) .*= 2

    for ℓ = 1:L
        Δℓ = Δ[ℓ]
        TwoLp1div2 = (2ℓ + 1) / 2

        for m = (-ℓ):ℓ
            m_ind = mod(m, 2 * L + 1) + 1
            Δm_ind = mod(m, 2 * ℓ + 1) + 1
            k = ℓ * (ℓ + 1) + m

            for m_ = 0:ℓ
                m__ind = mod(m_, 2 * L + 1) + 1
                Δm__ind = mod(m_, 2 * ℓ + 1) + 1
                for μ = -1:2:1
                    Δμ_ind = -1
                    μ_ind = mod(3 + μ, 3)
                    if μ == 1
                        Δμ_ind = 2
                    elseif μ == -1
                        Δμ_ind = 2 * ℓ + 1
                    end
                    u[k, μ_ind] +=
                        _complexunitpower(μ - m) *
                        Δℓ[Δm__ind, Δμ_ind] *
                        Δℓ[Δm__ind, Δm_ind] *
                        K[m__ind, m_ind, μ_ind] *
                        TwoLp1div2

                end

            end

        end

    end
    return u
end

function _θintegral!(u::AbstractArray{C}, v, dvec, weightvec, sinmθ, cosmθ, Δ, L) where {C}
    #irregularly distributed θ
    fill!(u, zero(C))
    for ℓ = 1:L
        Δℓ = Δ[ℓ]

        for m = (-ℓ):ℓ
            m_ind = mod(m, 2 * L + 1) + 1
            Δm_ind = mod(m, 2 * ℓ + 1) + 1
            k = ℓ * (ℓ + 1) + m
            for μ = -1:2:1
                Δμ_ind = -1
                μ_ind = mod(3 + μ, 3)
                if μ == 1
                    Δμ_ind = 2
                elseif μ == -1
                    Δμ_ind = 2 * ℓ + 1
                end
                mθind = 1
                ΔΔ = Δℓ[mθind, Δμ_ind] * Δℓ[mθind, Δm_ind]
                dvec .= ΔΔ
                trig = isodd(μ + m) ? sinmθ : cosmθ
                for mθ = 1:ℓ
                    mθind = mod(mθ, 2 * ℓ + 1) + 1
                    ΔΔ = Δℓ[mθind, Δμ_ind] * Δℓ[mθind, Δm_ind]
                    dvec .+= ΔΔ .* trig[mθ]
                end
                # dvec .*= (1.0im)^(μ-m)

                # u[k, μ_ind] += (2ℓ + 1) / 2 .* sum(dvec .* weightvec .* v[:, m_ind, μ_ind])
                dvec .*= _complexunitpower(μ - m) .* weightvec
                u[k, μ_ind] += (2ℓ + 1) ./ 2 .* (transpose(dvec) * view(v, :, m_ind, μ_ind))
            end

        end


    end
    return u
end


# using LinearAlgebra: det
# function u_to_β(u, αin,L)
function _receivecoeffs_2by2matrix(
    u::AbstractArray{C,2},
    αin::AbstractVector{C},
    L::Integer,
) where {C<:Complex}

    # βaut = zeros(C, 2 * L * (L + 2))
    βaut = Array{C}(undef, 2 * L * (L + 2))
    Amat = Matrix{C}(undef, 2, 2)
    uvectmp = Vector{C}(undef, 2)
    return _receivecoeffs_2by2matrix!(βaut, u, αin, L, Amat, uvectmp)
end
function _receivecoeffs_2by2matrix!(
    βaut::AbstractVector{C},
    u::AbstractArray{C,2},
    αin::AbstractVector{C},
    L::Integer,
    Amat,
    uvectmp,
) where {C<:Complex}
    for ℓ = 1:L
        Amat[1, 1] = αin[sℓm_to_j(1, ℓ, 1)]
        Amat[1, 2] = αin[sℓm_to_j(2, ℓ, 1)]
        Amat[2, 1] = αin[sℓm_to_j(1, ℓ, -1)]
        Amat[2, 2] = αin[sℓm_to_j(2, ℓ, -1)]

        for m = (-ℓ):ℓ
            βind = sℓm_to_j(1, ℓ, m)
            j = ℓ * (ℓ + 1) + m
            uvectmp[1], uvectmp[2] = u[j, 1], u[j, 2]
            # βaut[βind:(βind+1)] .= if abs(det(Amat)) > 1e-15
            βaut[βind:(βind+1)] .= if norm(Amat) > 1e-15
                SMatrix{2,2,C}(Amat) \ uvectmp
            else
                zeros(C, 2)
            end
        end
    end
    return βaut
end

# function αβ_to_u(α_inc,  β_aut,L;firstorder=true)
function _sum_product_αβ(
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    u = Vector{C}(undef, L * (L + 2), 2 * L + 1)
    return _sum_product_αβ!(u, α_inc, β_aut, L, firstorder = firstorder)
end

function _sum_product_αβ!(
    u,
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    ## |  u_μℓm = Σ_s α_sℓμ β_sℓm
    ## v
    # iterate over all ℓ m μ
    fill!(u, zero(C))
    for ℓ = 1:L
        for m = (-ℓ):ℓ
            k = ℓ * (ℓ + 1) + m
            # j = 2 * (ℓ * (ℓ + 1) + m - 1) + 1
            j = sℓm_to_j(1, ℓ, m)
            μrange = firstorder ? (-1:2:1) : (-ℓ:ℓ)

            for μ in μrange
                # μind=mod(3+μ,3)
                μind = firstorder ? mod(3 + μ, 3) : mod(μ, 2 * L + 1) + 1

                # j_ = 2 * (ℓ * (ℓ + 1) + μ - 1) + 1
                j_ = sℓm_to_j(1, ℓ, μ)
                u[k, μind] += β_aut[j] * α_inc[j_]
                u[k, μind] += β_aut[j+1] * α_inc[j_+1]
            end
        end
    end
    ## ^
    ## | u_μℓm  Σ_s α_sℓμ β_sℓm
    return u
end

# function αβ_to_u_ad(u,α_inc,L;firstorder=true)
function _sum_product_αβ_ad!(
    u,
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    ## |  u_μℓm = Σ_s α_sℓμ β_sℓm
    ## v 
    fill!(β_aut, zero(C))
    # iterate over all ℓ m μ
    for ℓ = 1:L
        for m = (-ℓ):ℓ

            μrange = -1:2:1
            if !firstorder
                μrange = (-ℓ):ℓ
            end
            for μ in μrange
                #   μind=mod(3+μ,3)
                μind = mod(3 + μ, 3)
                if !firstorder
                    μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                end
                k = ℓ * (ℓ + 1) + m
                j = 2 * (ℓ * (ℓ + 1) + m - 1) + 1
                j_ = 2 * (ℓ * (ℓ + 1) + μ - 1) + 1
                β_aut[j] += u[k, μind] * conj.(α_inc[j_])
                β_aut[j+1] += u[k, μind] * conj.(α_inc[j_+1])
            end
        end
    end
    ## ^
    ## | u_μℓm  Σ_s α_sℓμ β_sℓm
    return β_aut
end

# function u_to_v_(u,L;firstorder=true)
function _expandθmodes(
    u::AbstractArray{C,2},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    Nχ = 2
    Δ = Vector{Matrix{Float64}}(undef, L)
    for ℓ = 1:L
        Δ[ℓ] = ifftshift!(ifftshift!(_Δℓ_mμ(ℓ), 1), 2)
    end
    v_ = Vector{C}(undef, 2 * L + 1, 2 * L + 1, Nχ)

    return _expandθmodes!(v_, u, L, Δ; firstorder = firstorder)
end
function _expandθmodes!(
    v_,
    u::AbstractArray{C,2},
    L::Integer,
    Δ;
    firstorder = true,
) where {C<:Complex}
    ## |  v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    ## v

    fill!(v_, zero(C))
    # reset=true
    for ℓ = L:-1:1
        Δℓ = Δ[ℓ]
        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1
            Δmind = mod(m, 2 * ℓ + 1) + 1
            k = ℓ * (ℓ + 1) + m
            for mθ = (-ℓ):ℓ
                mθind = mod(mθ, 2 * L + 1) + 1 #mθ+L+1
                Δmθind = mod(mθ, 2 * ℓ + 1) + 1 #mθ+ℓ+1
                if firstorder
                    for μ = -1:2:1

                        if μ == 1
                            Δμind = 2
                            μind = 1
                        elseif μ == -1
                            Δμind = 2 * ℓ + 1
                            μind = 2
                        end

                        v_[mθind, mind, μind] +=
                            Δℓ[Δmθind, Δμind] *
                            Δℓ[Δmθind, Δmind] *
                            u[k, μind] *
                            _complexunitpower(m - μ)
                        # (1.0im)^(m - μ)


                    end
                else
                    for μ = (-ℓ):ℓ
                        μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                        Δμind = mod(μ, 2 * ℓ + 1) + 1#μ+ℓ+1
                        v_[mθind, mind, μind] +=
                            Δℓ[Δmθind, Δμind] *
                            Δℓ[Δmθind, Δmind] *
                            u[k, μind] *
                            _complexunitpower(m - μ)
                        # (1.0im)^(m - μ)

                    end
                end
            end
        end
    end
    ## ^
    ## | v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    return v_
end

function _complexunitpower(n::I) where {I<:Integer}
    pow = mod(n, 4)
    if pow == 0
        return Complex{I}(1, 0)
    elseif pow == 1
        return Complex{I}(0, 1)
    elseif pow == 2
        return Complex{I}(-1, 0)
    elseif pow == 3
        return Complex{I}(0, -1)
    end

end

function _expandθmodes_ad!(
    v_,
    u::AbstractArray{C,2},
    L::Integer,
    Δ;
    firstorder = true,
) where {C<:Complex}
    ## |  v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    ## v
    fill!(u, zero(C))
    for ℓ = 1:L
        Δℓ = Δ[ℓ]
        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1
            Δmind = mod(m, 2 * ℓ + 1) + 1
            k = ℓ * (ℓ + 1) + m
            for mθ = (-ℓ):ℓ
                mθind = mod(mθ, 2 * L + 1) + 1 #mθ+L+1
                Δmθind = mod(mθ, 2 * ℓ + 1) + 1 #mθ+ℓ+

                if firstorder
                    for μ = -1:2:1
                        Δμind = -1
                        μind = mod(3 + μ, 3) # μind=1 => μ=1;   μind=2 => μ=-1
                        if μ == 1
                            Δμind = 2
                        elseif μ == -1
                            Δμind = 2 * ℓ + 1
                        end
                        u[k, μind] +=
                            conj.(
                                Δℓ[Δmθind, Δμind] *
                                Δℓ[Δmθind, Δmind] *
                                _complexunitpower(m - μ)
                            ) * v_[mθind, mind, μind]

                    end
                else
                    for μ = (-ℓ):ℓ
                        μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                        Δμind = mod(μ, 2 * ℓ + 1) + 1#μ+ℓ+1
                        u[k, μind] +=
                            conj.(
                                Δℓ[Δmθind, Δμind] *
                                Δℓ[Δmθind, Δmind] *
                                _complexunitpower(m - μ)
                            ) * v_[mθind, mind, μind]

                    end
                end
            end
        end
    end
    ## ^
    ## | v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    return u
end


function _upsampleθ!(V, v_, L)
    view(V, 1:(L+1), :, :) .= view(v_, 1:(L+1), :, :)
    view(V, (size(V, 1)-L+1):size(V, 1), :, :) .= view(v_, (L+2):size(v_, 1), :, :)
    view(V, (L+2):(size(V, 1)-L), :, :) .= 0

    return V
end
function _upsampleθ_ad!(V, v_, L)
    view(v_, 1:(L+1), :, :) .= view(V, 1:(L+1), :, :)
    view(v_, (L+2):size(v_, 1), :, :) .= view(V, (size(V, 1)-L+1):size(V, 1), :, :)
    return v_
end
function _downsampleθ!(V, v_, Lθ, Jθ)
    V[1:(Lθ+1), :, :] = v_[1:(Lθ+1), :, :]
    if isodd(Jθ)
        V[(Lθ+2):end, :, :] = v_[(end-Lθ+1):end, :, :]
    else
        V[(Lθ+2):end, :, :] = v_[(end-Lθ):end, :, :]
    end
    return V
end
function _zeropaddingθ!(
    V::AbstractArray{C},
    v_::AbstractArray{C,3},
    L::Integer,
    Jθ::Integer,
    Lθ::Integer,
) where {C<:Complex}
    if (2 * L + 1) == Jθ
        V .= v_
        return V
    end
    if 2 * L + 1 < Jθ
        return _upsampleθ!(V, v_, L)

    elseif 2 * L + 1 > Jθ
        return _downsampleθ!(V, v_, Lθ, Jθ)
    end

end
function _zeropaddingθ(
    v_::AbstractArray{C,3},
    L::Integer,
    Jθ::Integer,
    Lθ::Integer;
    firstorder = true,
) where {C<:Complex}
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    V = Array{C}(undef, Jθ, 2 * L + 1, Nχ)
    return zeropaddingθ!(V, v_, L, Jθ, Lθ)
end
function _zeropaddingθ_ad!(
    V::AbstractArray{C},
    v_::AbstractArray{C,3},
    L::Integer,
    Jθ::Integer,
    Lθ::Integer,
) where {C<:Complex}
    if (2 * L + 1) == Jθ
        v_ .= V
        return v_
    end
    if 2 * L + 1 < Jθ
        v_ .= _upsampleθ_ad!(V, v_, L)
        return v_

    elseif 2 * L + 1 > Jθ
        v_ .= _downsampleθ_ad!(V, v_, Lθ, Jθ)
        v_
    end

end


function _upsampleϕ!(V, v, L, Jϕ)
    stopindex = size(v, 2)
    view(V, :, 1:(L+1), :) .= view(v, :, 1:(L+1), :)
    V[:, (Jϕ-L+1):end, :] .= view(v, :, (L+2):stopindex, :)
    V[:, (L+2):(end-L), :] .= 0
    return V
end
function _upsampleϕ_ad!(V, v, L, Jϕ)
    view(v, :, 1:(L+1), :) .= view(V, :, 1:(L+1), :)
    view(v, :, (L+2):size(v, 2), :) .= view(V, :, (Jϕ-L+1):size(V, 2), :)
    return v
end
function _downsampleϕ!(V, v, Jϕ)
    if isodd(Jϕ)
        Lϕ = Int((Jϕ - 1) / 2)
        V[:, 1:(Lϕ+1), :] .= v[:, 1:(Lϕ+1), :]
        V[:, (Lϕ+2):end, :] .= v[:, (end-Lϕ+1):end, :]
    else
        Lϕ = Int((Jϕ) / 2 - 1)
        V[:, 1:(Lϕ+1), :] .= v[:, 1:(Lϕ+1), :]
        V[:, (Lϕ+2):end, :] .= v[:, (end-Lϕ):end, :]
    end
    return V
end
function _zeropaddingϕ!(
    V::AbstractArray{C,3},
    v::AbstractArray{C,3},
    L::Integer,
    Jϕ::Integer,
) where {C<:Complex}
    if (2 * L + 1) == Jϕ
        view(V, 1:size(v, 1), 1:size(v, 2), :) .= v
        return V
    end

    if 2 * L + 1 < Jϕ
        V .= _upsampleϕ!(V, v, L, Jϕ)
        return V
    elseif 2 * L + 1 > Jϕ
        V .= _downsampleϕ!(V, v, Jϕ)
        return V
    end

end
function _zeropaddingϕ_ad!(
    V::AbstractArray{C,3},
    v::AbstractArray{C,3},
    L::Integer,
    Jϕ::Integer,
) where {C<:Complex}
    fill!(v, zero(C))
    if (2 * L + 1) == Jϕ
        v .= view(V, 1:size(v, 1), 1:size(v, 2), :)
        return v
    end

    if 2 * L + 1 < Jϕ
        v .= _upsampleϕ_ad!(V, v, L, Jϕ)
        return v
    elseif 2 * L + 1 > Jϕ
        v .= _downsampleϕ_ad!(V, v, Jϕ)
        return v
    end

end



function _χmodes_to_S12!(S12, w::AbstractArray{C,3}; firstorder = true) where {C<:Complex}
    if firstorder
        view(S12, :, :, 1) .= view(w, :, :, 1) .+ view(w, :, :, 2)
        view(S12, :, :, 2) .= view(w, :, :, 2) .- view(w, :, :, 1)
        view(S12, :, :, 2) .*= C(0, 1)
        return S12
    else
        Nθ, Jϕ, Jχ = size(w)
        oversamplingfactor = Jχ ÷ 4 + 1
        L = (Jχ - 1) ÷ 2
        v = zeros(C, Nθ, Jϕ, oversamplingfactor * 4)
        v = _zeropaddingχ!(v, w, L, oversamplingfactor * 4)
        S12 .= fft(v, 3)[:, :, [1, oversamplingfactor + 1]]
        return S12
    end
end
function _χmodes_to_S12_ad!(
    S12,
    w::AbstractArray{C,3};
    firstorder = true,
) where {C<:Complex}
    if firstorder
        view(w, :, :, 1) .= C(0, 1) .* view(S12, :, :, 2)
        view(w, :, :, 1) .+= view(S12, :, :, 1)
        view(w, :, :, 2) .= C(0, -1) .* view(S12, :, :, 2)
        view(w, :, :, 2) .+= view(S12, :, :, 1)
        return w
    else
        Nθ, Jϕ, Jχ = size(w)
        oversamplingfactor = Jχ ÷ 4 + 1
        L = (Jχ - 1) ÷ 2
        v = zeros(C, Nθ, Jϕ, oversamplingfactor * 4)
        view(v, :, :, [1, oversamplingfactor + 1]) .= S12
        v = bfft(v, 3)
        w .= _zeropaddingχ_ad!(v, w, L, oversamplingfactor * 4)
        return w
    end
end

function _inputdimensions(α_inc, β_aut, Jθ)
    # Determine dimensions of input data
    minlength = min(length(α_inc), length(β_aut))
    s = 0
    if isodd(minlength)
        s = 1
    else
        s = 2
    end
    # L is the maximum mode order
    L = Int(floor(sqrt((minlength - s) / 2 + 1)))

    Nθ = 0
    Lθ = 0
    if isodd(Jθ)
        Lθ = Int((Jθ - 1) / 2)
        Nθ = Lθ + 1
    else
        Lθ = Int(Jθ / 2 - 1)
        Nθ = Lθ + 2
    end
    return L, Nθ, Lθ
end


function _storage_fastspherical(
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    Jθ,
    Jϕ;
    firstorder = true,
) where {C}
    L, Nθ, Lθ = _inputdimensions(α_inc, β_aut, Jθ)
    # A1= L * (L + 2) * (2 * L + 1)
    # A2 = Jθ*(2L+1) *2
    # A3 = Nθ * Jϕ *2
    # maxA= maximum([A1, A2, A3])

    # B1= (2 * L + 1)*(2 * L + 1)* 2
    # B2= (L * (L + 2)) * (2 * L + 1) *2
    # B3 = Nθ * Jϕ *2
    # maxB= maximum([B1, B2, B3])

    # Atmp= zeros(C, maxA)
    # Btmp= zeros(C, maxB)

    Jθoversampled = Jθ
    oversamplingfactorθ = 1
    if 2 * L + 1 > Jθ
        oversamplingfactorθ = (2 * L + 1 ÷ Jθ) + 1
        Jθoversampled = oversamplingfactorθ * Jθ
    end

    Jϕoversampled = Jϕ
    oversamplingfactorϕ = 1
    if 2 * L + 1 > Jϕ
        oversamplingfactorϕ = (2 * L + 1 ÷ Jϕ) + 1
        Jϕoversampled = oversamplingfactorϕ * Jϕ
    end

    indθ = L * (L + 2)
    indϕ = 2 * L + 1
    # u= reshape(view(Atmp, 1: (indθ * indϕ)), indθ, indϕ)
    u = zeros(C, indθ, indϕ)

    Nχ = firstorder ? 2 : 2 * L + 1

    indθ = 2 * L + 1
    indϕ = 2 * L + 1
    # v__ = reshape(view(Btmp, 1: (2 * indθ * indϕ) ), indθ, indϕ, 2)
    v__ = zeros(C, indθ, indϕ, Nχ)

    indθ = Jθoversampled
    indϕ = 2 * L + 1
    # v_ = reshape(view(Atmp, 1: (2 * indθ * indϕ)), indθ, indϕ, 2)
    v_ = zeros(C, indθ, indϕ, Nχ)

    indθ = Nθ
    indϕ = Jϕoversampled
    # v = reshape(view(Btmp, 1: (2 * indθ * indϕ) ), indθ, indϕ, 2)
    v = zeros(C, indθ, indϕ, Nχ)
    # S21 = reshape(view(Atmp, 1: (2 * indθ * indϕ)), indθ, indϕ, 2)
    S21 = zeros(C, indθ, Jϕ, 2)

    Δ = Vector{Matrix{Float64}}(undef, L)
    for ℓ = 1:L
        Δ[ℓ] = Matrix{Float64}(undef, 2ℓ + 1, 2ℓ + 1)
        ifftshift!(Δ[ℓ], _Δℓ_mμ(ℓ))
    end

    fftplanθ! = plan_fft!(v_, 1)
    fftplanϕ! = plan_fft!(v, 2)


    return L,
    Nθ,
    Lθ,
    u,
    v__,
    v_,
    v,
    S21,
    Δ,
    fftplanθ!,
    fftplanϕ!,
    Jθoversampled,
    Jϕoversampled

end

function _expandirregularθ!(L, u, v_, cosmθ, sinmθ, Δ)
    for ℓ = 1:L

        Δℓ = Δ[ℓ]

        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1 #m+L+1
            Δmind = mod(m, 2 * ℓ + 1) + 1 #m+ℓ+1
            k = ℓ * (ℓ + 1) + m
            for μ = -1:2:1
                Δμind = -1
                μind = mod(3 + μ, 3) # μind=1 => μ=1;   μind=2 => μ=-1
                if μ == 1
                    Δμind = 2
                elseif μ == -1
                    Δμind = 2 * ℓ + 1
                end
                mθind = 1
                ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]
                imfac = (_complexunitpower(μ - m) * u[k, μind])

                # imfac = ((1.0im)^(μ - m) * u[k, μind])

                vview = view(v_, :, mind, μind)
                vview .+= imfac .* ΔΔ

                trig = isodd(μ + m) ? sinmθ : cosmθ
                for mθ = 1:ℓ
                    mθind = mod(mθ, 2 * ℓ + 1) + 1
                    ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]

                    vview .+= imfac .* ΔΔ .* trig[mθ]
                end


            end
        end
    end
    return v_
end
function _expandirregularθ_ad!(L, u, v_, cosmθ, sinmθ, Δ)
    fill!(u, zero(eltype(u)))
    for ℓ = 1:L

        Δℓ = Δ[ℓ]

        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1 #m+L+1
            Δmind = mod(m, 2 * ℓ + 1) + 1 #m+ℓ+1
            k = ℓ * (ℓ + 1) + m
            for μ = -1:2:1
                Δμind = -1
                μind = mod(3 + μ, 3) # μind=1 => μ=1;   μind=2 => μ=-1
                if μ == 1
                    Δμind = 2
                elseif μ == -1
                    Δμind = 2 * ℓ + 1
                end
                mθind = 1
                ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]
                imfac = (_complexunitpower(μ - m))

                # imfac = ((1.0im)^(μ - m) * u[k, μind])

                vview = view(v_, :, mind, μind)
                u[k, μind] += sum(conj.(imfac .* ΔΔ) .* vview)

                trig = isodd(μ + m) ? sinmθ : cosmθ
                for mθ = 1:ℓ
                    mθind = mod(mθ, 2 * ℓ + 1) + 1
                    ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]

                    u[k, μind] += sum(conj.(imfac .* ΔΔ .* trig[mθ]) .* vview)
                end


            end
        end
    end
    return u
end
function _storage_fastspherical_irregularθ(
    αin::AbstractVector{C},
    β_aut::AbstractVector{C},
    θvec,
    Jϕ,
) where {C}
    nthet = length(θvec)
    s = 0
    if isodd(length(αin))
        s = 1
    else
        s = 2
    end
    # L is the maximum mode order
    L, Nθ, Lθ = _inputdimensions(αin, β_aut, 2 * length(θvec))

    Jϕoversampled = Jϕ
    oversamplingfactorϕ = 1
    if 2 * L + 1 > Jϕ
        oversamplingfactorϕ = (2 * L + 1 ÷ Jϕ) + 1
        Jϕoversampled = oversamplingfactorϕ * Jϕ
    end


    indθ = L * (L + 2)
    indϕ = 2 * L + 1
    # u= reshape(view(Atmp, 1: (indθ * indϕ)), indθ, indϕ)
    u = zeros(C, indθ, indϕ)

    v_ = zeros(C, nthet, indϕ, 2)
    indθ = nthet
    indϕ = Jϕoversampled
    # v = reshape(view(Btmp, 1: (2 * indθ * indϕ) ), indθ, indϕ, 2)
    v = zeros(C, indθ, indϕ, 2)
    S21 = zeros(C, indθ, Jϕ, 2)


    cosmθ = [2 * cos.(mθ * θvec) for mθ = 1:L]
    sinmθ = [complex(0, 2) * sin.(mθ * θvec) for mθ = 1:L]
    # indices = [ifftshift(1:(2ℓ+1)) for ℓ = 1:L]
    Δ = Vector{Matrix{Float64}}(undef, L)
    for ℓ = 1:L
        Δ[ℓ] = Matrix{Float64}(undef, 2ℓ + 1, 2ℓ + 1)
        ifftshift!(Δ[ℓ], _Δℓ_mμ(ℓ))
    end
    fftplanϕ! = plan_fft!(v, 2)



    return L, u, v_, v, S21, cosmθ, sinmθ, Δ, fftplanϕ!, Jϕoversampled
end
function _zeropaddingχ!(
    W::AbstractArray{C,3},
    w::AbstractArray{C,3},
    L::Integer,
    Jχ::Integer,
) where {C<:Complex}
    if 2 * L + 1 == Jχ
        W .= w
        return W
    end
    if 2 * L + 1 > Jχ
        if isodd(Jχ)
            Lχ = Int((Jχ - 1) / 2)
            view(W, :, :, 1:(Lχ+1)) .= view(w, :, :, 1:(Lχ+1))
            view(W, :, :, ((Lχ+2):size(W, 3))) .=
                view(w, :, :, (size(w, 3)-Lχ+1):size(w, 3))
            return W
        else
            Lχ = Int((Jχ) / 2 - 1)
            view(W, :, :, 1:(Lχ+1)) .= view(w, :, :, 1:(Lχ+1))
            view(W, :, :, ((Lχ+2):size(W, 3))) .=
                view(w, :, :, ((size(w, 3)-Lχ):size(w, 3)))
            return W
        end
    elseif 2 * L + 1 < Jχ
        view(W, :, :, 1:(L+1)) .= view(w, :, :, 1:(L+1))
        view(W, :, :, (Jχ-L+1):size(W, 3)) .= view(w, :, :, (L+2):size(w, 3))
        return W
    end
    return W
end
function _zeropaddingχ(
    w::AbstractArray{C,3},
    L::Integer,
    Jχ::Integer,
    Jϕ::Integer,
    Nθ::Integer,
) where {C<:Complex}
    W = zeros(C, Nθ, Jϕ, Jχ)
    return _zeropaddingχ!(W, w, L, Jχ)
end
function _zeropaddingχ_ad!(
    W::AbstractArray{C,3},
    w::AbstractArray{C,3},
    L::Integer,
    Jχ::Integer,
) where {C<:Complex}
    if 2 * L + 1 == Jχ
        w .= W
        return w
    end
    if 2 * L + 1 > Jχ
        if isodd(Jχ)
            Lχ = Int((Jχ - 1) / 2)
            view(w, :, :, 1:(Lχ+1)) .= view(W, :, :, 1:(Lχ+1))
            view(w, :, :, (size(w, 3)-Lχ+1):size(w, 3)) .=
                view(W, :, :, ((Lχ+2):size(W, 3)))
            return w
        else
            Lχ = Int((Jχ) / 2 - 1)
            view(w, :, :, 1:(Lχ+1)) .= view(W, :, :, 1:(Lχ+1))
            view(w, :, :, ((size(w, 3)-Lχ):size(w, 3))) .=
                view(W, :, :, ((Lχ+2):size(W, 3)))
            return w
        end
    elseif 2 * L + 1 < Jχ
        view(w, :, :, 1:(L+1)) .= view(W, :, :, 1:(L+1))
        view(w, :, :, (L+2):size(w, 3)) .= view(W, :, :, (Jχ-L+1):size(W, 3))
        return w
    end
    return w
end



function _φ_integral!(v, w, L, ifft_planϕ!)
    #irregularly distributed θ
    mul!(w, ifft_planϕ!, w)
    view(v, :, 1:L+1, :) .= view(w, :, 1:(L+1), :)  # remove redundant sample of ϕ-harmonics
    view(v, :, L+2:size(v, 2), :) .= view(w, :, (size(w, 2)-L+1):size(w, 2), :) # remove redundant sample of ϕ-harmonics
    return v
end

function _storage_fastsphericalinverse(
    S12::AbstractArray{C,3},
    θvec::AbstractArray{F,1},
) where {C,F}
    #irregularly distributed θ 

    a, b, _ = size(S12)

    L = minimum([a - 1, (b - 1) ÷ 2])
    w = Array{C}(undef, a, b, 2)
    u = zeros(C, L * (L + 2), 2)
    v = Array{C}(undef, a, 2L + 1, 2)
    nthet = length(θvec)
    cosmθ = [2 * cos.(mθ * θvec) for mθ = 1:L]
    sinmθ = [complex(0, 2) * sin.(mθ * θvec) for mθ = 1:L]
    dvec = Array{C}(undef, nthet)
    Δ = Vector{Matrix{Float64}}(undef, L)
    for ℓ = 1:L
        Δ[ℓ] = Matrix{Float64}(undef, 2ℓ + 1, 2ℓ + 1)
        ifftshift!(Δ[ℓ], _Δℓ_mμ(ℓ))
    end

    βaut = Vector{C}(undef, 2 * L * (L + 2))
    αaut = Vector{C}(undef, 2 * L * (L + 2))
    Amat = Matrix{C}(undef, 2, 2)
    uvectmp = Vector{C}(undef, 2)
    ifft_planϕ! = plan_ifft!(w, 2)

    return w, u, v, L, ifft_planϕ!, cosmθ, sinmθ, dvec, βaut, αaut, Amat, uvectmp, Δ
end

function _storage_fastsphericalinverse(
    S12::AbstractArray{C,3},
    Jθ::Integer,
    Jϕ::Integer,
) where {C}
    Lϕ = (Jϕ - 1) ÷ 2
    Lθ = (Jθ - 1) ÷ 2
    ## | Determine size of input data
    ## v   
    a, b, c = size(S12)
    Nθexpected = isodd(Jθ) ? Lθ + 1 : Lθ + 2
    if Nθexpected != (Jθ - Lθ)
        println("Ohoh")
    end
    a == Nθexpected ? () :
    throw(DimensionMismatch("First dimension of S12 is not consistent with Jθ."))
    b == Jϕ ? () :
    throw(DimensionMismatch("Second dimension of S12 is not consistent with Jϕ."))
    c == 2 ? () : throw(DimensionMismatch("Third dimension of S12 must be 2."))
    ## ^
    ## | Determine size of input data
    L = minimum([Lθ, Lϕ])
    # Lmax= maximum([Lθ, Lϕ])
    # S12padded = Array{C}(undef, L + 1, 2 * L + 1, 2)
    v = zeros(C, Jθ, Jϕ, 2)
    vview = view(v, 1:(Jθ-Lθ), :, :)
    vview2 = view(v, :, 1:2L+1, :)

    # vview = view(v, 1:(Lθ+1), :, :)
    ifft_planϕ = plan_ifft!(vview, 2)
    ifft_planθ = plan_ifft!(vview2, 1)
    u = zeros(C, L * (L + 2), 2)
    P = fft(_Πfun.(ifftshift((-2*L):(2*L))))
    K = zeros(C, 4 * L + 1, 2 * L + 1, 2)
    βaut = Array{C}(undef, 2 * L * (L + 2))
    # αin_tmp = zeros(C, 2 * L * (L + 2))
    αaut = Array{C}(undef, 2 * L * (L + 2))
    Amat = Matrix{C}(undef, 2, 2)
    uvectmp = Vector{C}(undef, 2)

    Δ = Vector{Matrix{Float64}}(undef, L)
    for ℓ = 1:L
        Δ[ℓ] = Matrix{Float64}(undef, 2ℓ + 1, 2ℓ + 1)
        ifftshift!(Δ[ℓ], _Δℓ_mμ(ℓ))
    end
    return L,
    v,
    vview,
    vview2,
    ifft_planϕ,
    ifft_planθ,
    u,
    P,
    K,
    βaut,
    αaut,
    Amat,
    uvectmp,
    Δ

end

####################################################################################################
#                         ^
#                         |
####################### auxilliary functions #######################################################


####################### fastsphericalforward #######################################################
#                         |
#                         v 
####################################################################################################
"""
    fastsphericalforward(α_inc::AbstractVector{C}, β_aut::AbstractVector{C}, ...)


Return S12 data with grid points on sphere from (receiving) AUT spherical mode coefficients 'β_aut' and  incident probe field coefficients 'α_inc'.
"""
function fastsphericalforward(
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    Jθ::Integer,
    Jϕ::Integer;
    firstorder = true,
) where {C<:Complex}
    # First-order or arbitrary order probe, regular sampling in θ

    L, Nθ, Lθ, u, v__, v_, v, S21, Δ, fftplanθ!, fftplanϕ!, Jθoversampled, Jϕoversampled =
        _storage_fastspherical(α_inc, β_aut, Jθ, Jϕ, firstorder = firstorder)

    return fastsphericalforward!(
        α_inc,
        β_aut,
        Jθ,
        Jϕ,
        Jθoversampled,
        Jϕoversampled,
        L,
        Nθ,
        Lθ,
        u,
        v__,
        v_,
        v,
        S21,
        Δ,
        fftplanθ!,
        fftplanϕ!,
        firstorder = firstorder,
    )
end
function fastsphericalforward(
    αin::AbstractVector{C},
    β_aut::AbstractVector{C},
    θvec::AbstractVector{F},
    Jϕ::Integer,
) where {C<:Complex,F<:Number}
    # First-order probe, irregular sampling in θ

    L, u, v_, v, S21, cosmθ, sinmθ, Δ, fftplanϕ!, Jϕoversampled =
        _storage_fastspherical_irregularθ(αin, β_aut, θvec, Jϕ)
    Nθ = length(θvec)
    return fastsphericalforward!(
        αin,
        β_aut,
        Nθ,
        Jϕ,
        Jϕoversampled,
        L,
        u,
        v_,
        v,
        S21,
        cosmθ,
        sinmθ,
        Δ,
        fftplanϕ!,
    )
end

"""
    fastsphericalforward!(α_inc::AbstractVector{C}, β_aut::AbstractVector{C}, ...)


Return S12 data with grid points on sphere from (receiving) AUT spherical mode coefficients 'β_aut' and  incident probe field coefficients 'α_inc'.
In-place version of `fastsphericalforward`.
"""
function fastsphericalforward!(
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    Jθ::Integer,
    Jϕ::Integer,
    Jθoversampled,
    Jϕoversampled,
    L,
    Nθ,
    Lθ,
    u,
    v__,
    v_,
    v,
    S21,
    Δ,
    fftplanθ!,
    fftplanϕ!;
    firstorder = true,
) where {C<:Complex}
    # First-order or arbitrary order probe, regular sampling in θ
    u .= _sum_product_αβ!(u, α_inc, β_aut, L; firstorder = firstorder)
    v__ .= _expandθmodes!(v__, u, L, Δ; firstorder = firstorder)


    # for correct downsample, we must first upsample to integer multiple of Jθ
    oversamplingfactorθ = Jθoversampled ÷ Jθ
    v_ .= _zeropaddingθ!(v_, v__, L, Jθoversampled, Lθ)
    mul!(v_, fftplanθ!, v_)

    # for correct downsample, we must first upsample to integer multiple of Jϕ
    oversamplingfactorϕ = Jϕoversampled ÷ Jϕ
    v .= _zeropaddingϕ!(
        v,
        view(v_, 1:oversamplingfactorθ:oversamplingfactorθ*Nθ, :, :),
        L,
        Jϕoversampled,
    )


    mul!(v, fftplanϕ!, v)

    return _χmodes_to_S12!(
        S21,
        view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :),
        firstorder = firstorder,
    )
end
function fastsphericalforward!(
    αin::AbstractVector{C},
    β_aut::AbstractVector{C},
    Nθ,
    Jϕ::Integer,
    Jϕoversampled,
    L,
    u,
    v_,
    v,
    S21,
    cosmθ,
    sinmθ,
    Δ,
    fftplanϕ!,
) where {C}
    # First-order probe, irregular sampling in θ

    u .= _sum_product_αβ!(u, αin, β_aut, L; firstorder = true)
    v_ .= _expandirregularθ!(L, u, v_, cosmθ, sinmθ, Δ)
    v .= _zeropaddingϕ!(v, v_, L, Jϕoversampled)
    mul!(v, fftplanϕ!, v)

    # for correct downsample, we must first upsample to integer multiple of Jϕ
    oversamplingfactorϕ = Jϕoversampled ÷ Jϕ

    return _χmodes_to_S12!(S21, view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :))
end



####################################################################################################
#                         ^
#                         |
####################### fastsphericalforward #######################################################

####################### fastsphericalforward adjoint ###############################################
#                         |
#                         v 
####################################################################################################


"""
    fastsphericalforward_ad!(
    αin::AbstractVector{C},
    β_aut::AbstractVector{C}, ...)
  
  
Return adjoint to `fastsphericalforward!(α_inc,  β_aut,  ...`.
"""
function fastsphericalforward_ad!(
    α_inc::AbstractVector{C},
    β_aut::AbstractVector{C},
    Jθ::Integer,
    Jϕ::Integer,
    Jθoversampled,
    Jϕoversampled,
    L,
    Nθ,
    Lθ,
    u,
    v__,
    v_,
    v,
    S21,
    Δ,
    fftplanθ!,
    fftplanϕ!;
    firstorder = true,
) where {C<:Complex}
    # First-order or arbitrary order probe, regular sampling in θ
    oversamplingfactorθ = Jθoversampled ÷ Jθ
    oversamplingfactorϕ = Jϕoversampled ÷ Jϕ
    fill!(v, zero(C))
    view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :) .= _χmodes_to_S12_ad!(
        S21,
        view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :),
        firstorder = firstorder,
    )

    bfftplanϕ! = (Jϕoversampled * inv(fftplanϕ!))
    mul!(v, bfftplanϕ!, v)
    fill!(v_, zero(C))
    view(v_, 1:oversamplingfactorθ:oversamplingfactorθ*Nθ, :, :) .= _zeropaddingϕ_ad!(
        v,
        view(v_, 1:oversamplingfactorθ:oversamplingfactorθ*Nθ, :, :),
        L,
        Jϕoversampled,
    )
    fftplanθ! = (Jθoversampled * inv(fftplanθ!))
    mul!(v_, fftplanθ!, v_)
    v__ .= _zeropaddingθ_ad!(v_, v__, L, Jθoversampled, Lθ)
    u .= _expandθmodes_ad!(v__, u, L, Δ; firstorder = firstorder)
    β_aut .= _sum_product_αβ_ad!(u, α_inc, β_aut, L; firstorder = firstorder)

    return β_aut
end
function fastsphericalforward_ad!(
    αin::AbstractVector{C},
    β_aut::AbstractVector{C},
    Nθ,
    Jϕ::Integer,
    Jϕoversampled,
    L,
    u,
    v_,
    v,
    S21,
    cosmθ,
    sinmθ,
    Δ,
    fftplanϕ!,
) where {C}
    # First-order probe, irregular sampling in θ
    oversamplingfactorϕ = Jϕoversampled ÷ Jϕ
    view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :) .= _χmodes_to_S12_ad!(
        S21,
        view(v, :, 1:oversamplingfactorϕ:Jϕoversampled, :),
        firstorder = true,
    )

    bfftplanϕ! = (UniformScaling(size(v, 2)) * inv(fftplanϕ!))
    mul!(v, bfftplanϕ!, v)
    v_ .= _zeropaddingϕ_ad!(v, v_, L, Jϕoversampled)
    u .= _expandirregularθ_ad!(L, u, v_, cosmθ, sinmθ, Δ)

    β_aut .= _sum_product_αβ_ad!(u, αin, β_aut, L; firstorder = true)

    return β_aut
end

####################################################################################################
#                         ^
#                         |
####################### fastsphericalforward adjoint ###############################################

####################### fastsphericalinverse  ######################################################
#                         |
#                         v 
####################################################################################################


"""
    fastsphericalinverse(S12, αin, ...)


Return (receiving) spherical mode coefficients from measured S12 data with grid points on sphere.
"""
function fastsphericalinverse(
    S12::AbstractArray{C,3},
    αin::AbstractArray{T,1},
    Jθ::Integer,
    Jϕ::Integer,
) where {C<:Complex} where {T<:Complex}


    L, v, vview, vview2, ifft_planϕ, ifft_planθ, u, P, K, βaut, αaut, Amat, uvectmp, Δ =
        _storage_fastsphericalinverse(S12, Jθ, Jϕ)
    return fastsphericalinverse!(
        S12,
        αin,
        Jθ,
        Jϕ,
        L,
        v,
        vview,
        vview2,
        ifft_planϕ,
        ifft_planθ,
        u,
        P,
        K,
        βaut,
        αaut,
        Amat,
        uvectmp,
        Δ,
    )
end
function fastsphericalinverse(
    S12::AbstractArray{C,3},
    αin::AbstractArray{C,1},
    θvec::AbstractArray{F,1},
    weightvec::AbstractArray{K,1},
) where {C<:Complex} where {F<:Real} where {K<:Number}

    w, u, v, L, ifft_planϕ!, cosmθ, sinmθ, dvec, βaut, αaut, Amat, uvectmp, Δ =
        _storage_fastsphericalinverse(S12, θvec)

    return fastsphericalinverse!(
        S12,
        αin,
        weightvec,
        L,
        w,
        u,
        v,
        ifft_planϕ!,
        cosmθ,
        sinmθ,
        dvec,
        βaut,
        αaut,
        Amat,
        uvectmp,
        Δ,
    )
end


"""
    fastsphericalinverse!(S12, αin, [θvec, weightvec])


Return spherical mode coefficients from measured S12 data with grid points on sphere.
In-place version of `fastsphericalinverse`.
"""
function fastsphericalinverse!(
    S12::AbstractArray{C,3},
    αin::AbstractArray{T,1},
    Jθ::Integer,
    Jϕ::Integer,
    L,
    v,
    vview,
    vview2,
    ifft_planϕ,
    ifft_planθ,
    u,
    P,
    K,
    βaut,
    αaut,
    Amat,
    uvectmp,
    Δ,
) where {C<:Complex} where {T<:Complex}

    vview .= _χ_integral!(vview, S12)
    v .= _φ_integral!(v, vview, L, ifft_planϕ, Jθ, Jϕ)
    u .= _θintegral!(u, vview2, L, P, K, Δ, ifft_planθ)
    return βtoα!(αaut, _receivecoeffs_2by2matrix!(βaut, u, αin, L, Amat, uvectmp))
end
function fastsphericalinverse!(
    S12::AbstractArray{C,3},
    αin::AbstractArray{C,1},
    weightvec::AbstractArray{K,1},
    L,
    w,
    u,
    v,
    ifft_planϕ!,
    cosmθ,
    sinmθ,
    dvec,
    βaut,
    αaut,
    Amat,
    uvectmp,
    Δ,
) where {C<:Complex} where {K<:Number}
    #irregularly distributed θ
    w .= _χ_integral!(w, S12)
    v .= _φ_integral!(v, w, L, ifft_planϕ!)
    u .= _θintegral!(u, v, dvec, weightvec, sinmθ, cosmθ, Δ, L)
    return βtoα!(αaut, _receivecoeffs_2by2matrix!(βaut, u, αin, L, Amat, uvectmp))
end



####################################################################################################
#                         ^
#                         |
####################### fastsphericalinverse  ######################################################
