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

function _Convfun(A, B, L)
    length(A) == 4 * L + 1 ? () : throw(DimensionMismatch("length(A) must be 4L+1"))
    length(B) == 2 * L + 1 ? () : throw(DimensionMismatch("length(B) must be 2L+1"))
    typevar = promote_type(eltype(A), eltype(B))
    return ifft(fft(A) .* fft([B[1:(L+1)]; zeros(typevar, 2 * L); B[(L+2):(2*L+1)]]))

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
function _χandφ_integral(S12::Array{C,3}, L::Integer) where {C<:Complex}
    v = zeros(C, 2 * L + 1, 2 * L + 1, 2)
    v[1:(L+1), 1:(2*L+1), 1] .= view(S12, :, :, 1) - C(0, 1) .* view(S12, :, :, 2) # w_1(ϑ, ϕ)
    v[1:(L+1), 1:(2*L+1), 2] .= view(S12, :, :, 1) + C(0, 1) .* view(S12, :, :, 2) # w_-1 (ϑ, ϕ)
    v[1:(L+1), 1:(2*L+1), :] .*= 0.5
    v[1:(L+1), :, :] .= ifft!(view(v, 1:(L+1), :, :), 2)

    v[(L+2):end, :, :] .= view(v, (L+1):-1:2, :, :)
    μval = fftfreq(2 * L + 1, 2 * L + 1)
    for k = 1:(2*L+1)
        signfac = (-1)^(μval[k] + 1)
        v[(L+2):end, k, :] .*= signfac
    end
    return v
end

# function v_to_u(v,L)
function _θintegral(v::Array{C,3}, L::Integer) where {C<:Number}

    P = _Πfun.(ifftshift((-2*L):(2*L)))
    K = zeros(C, 4 * L + 1, 2 * L + 1, 2)

    for m_ind = 1:(2*L+1)
        K[:, m_ind, 1] = _Convfun(P, v[:, m_ind, 1], L)
        K[:, m_ind, 2] = _Convfun(P, v[:, m_ind, 2], L)
    end

    K[2:end, :, :] *= 2

    u = zeros(C, L * (L + 2), 2)
    for ℓ = 1:L
        Δℓ = _Δℓ_mμ(ℓ)
        Δℓ = ifftshift(Δℓ, 1)
        Δℓ = ifftshift(Δℓ, 2)

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
                        (1.0im)^(μ - m) *
                        Δℓ[Δm__ind, Δμ_ind] *
                        Δℓ[Δm__ind, Δm_ind] *
                        K[m__ind, m_ind, μ_ind] *
                        (2ℓ + 1) / 2
                end

            end

        end

    end
    return u
end

using LinearAlgebra: det
# function u_to_β(u, αin,L)
function _receivecoeffs_2by2matrix(
    u::Array{C,2},
    αin::Vector{C},
    L::Integer,
) where {C<:Complex}

    # βaut = zeros(C, 2 * L * (L + 2))
    βaut = Array{C}(undef, 2 * L * (L + 2))
    for ℓ = 1:L
        Amat = [
            αin[sℓm_to_j(1, ℓ, 1)] αin[sℓm_to_j(2, ℓ, 1)]
            αin[sℓm_to_j(1, ℓ, -1)] αin[sℓm_to_j(2, ℓ, -1)]
        ]
        for m = (-ℓ):ℓ
            βind = sℓm_to_j(1, ℓ, m)
            j = ℓ * (ℓ + 1) + m
            βaut[βind:(βind+1)] .= if abs(det(Amat)) > 1e-15
                Amat \ [u[j, 1]; u[j, 2]]
            else
                zero(C)
            end
        end
    end
    return βaut
end

# function αβ_to_u(α_inc,  β_aut,L;firstorder=true)
function _sum_product_αβ(
    α_inc::Vector{C},
    β_aut::Vector{C},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    ## |  u_μℓm = Σ_s α_sℓμ β_sℓm
    ## v

    u = zeros(C, L * (L + 2), 2 * L + 1)
    # iterate over all ℓ m μ
    for ℓ = 1:L
        for m = (-ℓ):ℓ
            μrange = -1:2:1
            if !firstorder
                μrange = (-ℓ):ℓ
            end

            for μ in μrange
                # μind=mod(3+μ,3)
                μind = mod(3 + μ, 3)
                if !firstorder
                    μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                end
                k = ℓ * (ℓ + 1) + m
                j = 2 * (ℓ * (ℓ + 1) + m - 1) + 1
                j_ = 2 * (ℓ * (ℓ + 1) + μ - 1) + 1
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
function _sum_product_αβ_ad(
    u::Array{C,2},
    α_inc::Vector{C},
    L::Integer;
    firstorder = true,
) where {C<:Complex}
    ## |  u_μℓm = Σ_s α_sℓμ β_sℓm
    ## v 
    β_aut = zeros(C, size(α_inc))
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
                β_aut[j] += u[k, μind] * conj(α_inc[j_])
                β_aut[j+1] += u[k, μind] * conj(α_inc[j_+1])
            end
        end
    end
    ## ^
    ## | u_μℓm  Σ_s α_sℓμ β_sℓm
    return β_aut
end

# function u_to_v_(u,L;firstorder=true)
function _expandθmodes(u::Array{C,2}, L::Integer; firstorder = true) where {C<:Complex}
    ## |  v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    ## v
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    v_ = zeros(C, 2 * L + 1, 2 * L + 1, Nχ)
    for ℓ = 1:L
        Δℓ = _Δℓ_mμ(ℓ) # Wigner-d-matrices from WignerD.jl
        Δℓ = ifftshift(Δℓ, 1)
        Δℓ = ifftshift(Δℓ, 2)
        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1
            Δmind = mod(m, 2 * ℓ + 1) + 1
            k = ℓ * (ℓ + 1) + m
            for mθ = (-ℓ):ℓ
                mθind = mod(mθ, 2 * L + 1) + 1 #mθ+L+1
                Δmθind = mod(mθ, 2 * ℓ + 1) + 1 #mθ+ℓ+1
                if firstorder
                    for μ = -1:2:1
                        Δμind = -1
                        μind = mod(3 + μ, 3) # μind=1 => μ=1;   μind=2 => μ=-1
                        if μ == 1
                            Δμind = 2
                        elseif μ == -1
                            Δμind = 2 * ℓ + 1
                        end
                        v_[mθind, mind, μind] +=
                            Δℓ[Δmθind, Δμind] *
                            Δℓ[Δmθind, Δmind] *
                            u[k, μind] *
                            (1.0im)^(m - μ)

                    end
                else
                    for μ = (-ℓ):ℓ
                        μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                        Δμind = mod(μ, 2 * ℓ + 1) + 1#μ+ℓ+1
                        v_[mθind, mind, μind] +=
                            Δℓ[Δmθind, Δμind] *
                            Δℓ[Δmθind, Δmind] *
                            u[k, μind] *
                            (1.0im)^(m - μ)

                    end
                end
            end
        end
    end
    ## ^
    ## | v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    return v_
end
function _expandθmodes_ad(v_::Array{C,3}, L::Integer; firstorder = true) where {C<:Complex}
    ## |  v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    ## v
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    u = zeros(C, L * (L + 2), Nχ)
    for ℓ = 1:L
        Δℓ = _Δℓ_mμ(ℓ) # Wigner-d-matrices from WignerD.jl
        Δℓ = ifftshift(Δℓ, 1)
        Δℓ = ifftshift(Δℓ, 2)
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
                            conj(Δℓ[Δmθind, Δμind] * Δℓ[Δmθind, Δmind] * (1.0im)^(m - μ)) *
                            v_[mθind, mind, μind]

                    end
                else
                    for μ = (-ℓ):ℓ
                        μind = mod(μ, 2 * L + 1) + 1 #μ+L+1
                        Δμind = mod(μ, 2 * ℓ + 1) + 1#μ+ℓ+1
                        u[k, μind] +=
                            conj(Δℓ[Δmθind, Δμind] * Δℓ[Δmθind, Δmind] * (1.0im)^(m - μ)) *
                            v_[mθind, mind, μind]

                    end
                end
            end
        end
    end
    ## ^
    ## | v_{μ,m,mθ} = j^{m-μ} ∑_ℓ Δ^{ℓ}_mθ,μ Δ^{ℓ}_mθ,m u_μℓm
    return u
end



function _zeropaddingθ!(
    v_::Array{C,3},
    L::Integer,
    Jθ::Integer,
    Lθ::Integer;
    firstorder = true,
) where {C<:Complex}
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    if 2 * L + 1 < Jθ
        V = Array{C}(undef, Jθ, 2 * L + 1, Nχ)
        V[1:(L+1), :, :] = v_[1:(L+1), :, :]
        V[(end-L+1):end, :, :] = v_[(L+2):end, :, :]
        V[(L+2):(end-L), :, :] .= 0
        return V
        # V=0 # release memory

    elseif 2 * L + 1 > Jθ
        V = Array{C}(undef, Jθ, 2 * L + 1, Nχ)
        V[1:(Lθ+1), :, :] = v_[1:(Lθ+1), :, :]
        if isodd(Jθ)
            V[(Lθ+2):end, :, :] = v_[(end-Lθ+1):end, :, :]
        else
            V[(Lθ+2):end, :, :] = v_[(end-Lθ):end, :, :]
        end
        return V
        # V=0 # release memory
    end
    return v_
end
function _zeropaddingθ_ad(
    V_::Array{C,3},
    L::Integer,
    Jθ::Integer,
    Lθ::Integer;
    firstorder = true,
) where {C<:Complex}
    #remove zero-padding in θ
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    v_ = zeros(C, 2 * L + 1, 2 * L + 1, Nχ)
    if 2 * L + 1 < Jθ
        v_[1:(L+1), :, :] = V_[1:(L+1), :, :]
        v_[(L+2):end, :, :] = V_[(end-L+1):end, :, :]
        # V=0 # release memory

    elseif 2 * L + 1 > Jθ
        v_[1:(Lθ+1), :, :] = V_[1:(Lθ+1), :, :]
        if isodd(Jθ)
            v_[(end-Lθ+1):end, :, :] = V_[(Lθ+2):end, :, :]
        else
            v_[(end-Lθ):end, :, :] = V_[(Lθ+2):end, :, :]
        end
        # V=0 # release memory
    end
    return v_
end

function _zeropaddingϕ!(
    v::Array{C,3},
    L::Integer,
    Jϕ::Integer,
    Nθ::Integer;
    firstorder = true,
) where {C<:Complex}
    Nχ = 2
    if !firstorder
        Nχ = size(v, 3)
    end
    if 2 * L + 1 < Jϕ
        V = Array{C}(undef, Nθ, Jϕ, Nχ)
        V[:, 1:(L+1), :] = v[:, 1:(L+1), :]
        V[:, (Jϕ-L+1):end, :] = v[:, (L+2):end, :]
        V[:, (L+2):(Jϕ-L), :] .= 0
        return V
        # v = V
        # V=0 # release memory

    elseif 2 * L + 1 > Jϕ
        V = Array{C}(undef, Nθ, Jϕ, Nχ)
        if isodd(Jϕ)
            Lϕ = Int((Jϕ - 1) / 2)
            V[:, 1:(Lϕ+1), :] = v[:, 1:(Lϕ+1), :]
            V[:, (Lϕ+2):end, :] = v[:, (end-Lϕ+1):end, :]
        else
            Lϕ = Int((Jϕ) / 2 - 1)
            V[:, 1:(Lϕ+1), :] = v[:, 1:(Lϕ+1), :]
            V[:, (Lϕ+2):end, :] = v[:, (end-Lϕ):end, :]
        end
        return V
        # v = V
        # V=0 # release memory
    end
    return v
end
function _zeropaddingϕ_ad(
    V::Array{C,3},
    L::Integer,
    Jϕ::Integer,
    Nθ::Integer;
    firstorder = true,
) where {C<:Complex}
    Nχ = 2
    if !firstorder
        Nχ = 2 * L + 1
    end
    v = zeros(C, Jθ, 2 * L + 1, Nχ)
    if 2 * L + 1 < Jϕ
        v[1:Nθ, 1:(L+1), :] = V[:, 1:(L+1), :]
        v[1:Nθ, (L+2):end, :] = V[:, (Jϕ-L+1):end, :]
        # V=0 # release memory

    elseif 2 * L + 1 > Jϕ
        if isodd(Jϕ)
            Lϕ = Int((Jϕ - 1) / 2)
            v[1:Nθ, 1:(Lϕ+1), :] = V[:, 1:(Lϕ+1), :]
            v[1:Nθ, (end-Lϕ+1):end, :] = V[:, (Lϕ+2):end, :]
        else
            Lϕ = Int((Jϕ) / 2 - 1)
            v[1:Nθ, 1:(Lϕ+1), :] = V[:, 1:(Lϕ+1), :]
            v[1:Nθ, (end-Lϕ):end, :] = V[:, (Lϕ+2):end, :]
        end
        # V=0 # release memory
    else
        v[1:Nθ, :, :] = V
    end
    return v
end

function _zeropaddingχ!(
    w::Array{C,3},
    L::Integer,
    Jχ::Integer,
    Jϕ::Integer,
    Nθ::Integer,
) where {C<:Complex}
    if 2 * L + 1 > Jχ
        W = zeros(C, Nθ, Jϕ, Jχ)
        if isodd(Jχ)
            Lχ = Int((Jχ - 1) / 2)
            W[:, :, 1:(Lχ+1)] = w[:, :, 1:(Lχ+1)]
            W[:, :, (Lχ+2):end] = w[:, :, (end-Lχ+1):end]
        else
            Lχ = Int((Jχ) / 2 - 1)
            W[:, :, 1:(Lχ+1)] = w[:, :, 1:(Lχ+1)]
            W[:, :, (Lχ+2):end] = w[:, :, (end-Lχ):end]
        end
        w = W
    elseif 2 * L + 1 < Jχ
        W = zeros(C, Nθ, Jϕ, Jχ)
        W[:, :, 1:(L+1)] = w[:, :, 1:(L+1)]
        W[:, :, (Jχ-L+1):end] = w[:, :, (L+2):end]
        w = W
    end
    return w
end
function _zeropaddingχ_ad(
    W::Array{C,3},
    L::Integer,
    Jχ::Integer,
    Jϕ::Integer,
    Nθ::Integer,
) where {C<:Complex}
    w = zeros(C, Nθ, Jϕ, 2 * L + 1)
    if 2 * L + 1 > Jχ
        if isodd(Jχ)
            Lχ = Int((Jχ - 1) / 2)
            w[:, :, 1:(Lχ+1)] = W[:, :, 1:(Lχ+1)]
            w[:, :, (end-Lχ+1):end] = W[:, :, (Lχ+2):end]
        else
            Lχ = Int((Jχ) / 2 - 1)
            w[:, :, 1:(Lχ+1)] = W[:, :, 1:(Lχ+1)]
            w[:, :, (end-Lχ):end] = W[:, :, (Lχ+2):end]
        end
    elseif 2 * L + 1 < Jχ
        w[:, :, 1:(L+1)] = W[:, :, 1:(L+1)]
        w[:, :, (L+2):end] = W[:, :, (Jχ-L+1):end]
    else
        w = W
    end
    return w
end

# function w_to_S12_Firstorder(w)
function _χmodes_to_S12_Firstorder(w::Array{C,3}) where {C<:Complex}
    S12 = Array{C}(undef, size(w))
    S12[:, :, 1] .= w[:, :, 1] .+ w[:, :, 2]
    S12[:, :, 2] .= C(0, 1) .* (w[:, :, 1] .- w[:, :, 2])
    return S12
end

function _inputdimensions(α_inc, β_aut, Jθ)
    # Determine dimensions of input data
    s = 0
    if isodd(min(length(α_inc), length(β_aut)))
        s = 1
    else
        s = 2
    end
    # L is the maximum mode order
    L = Int(floor(sqrt((min(length(α_inc), length(β_aut)) - s) / 2 + 1)))

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

#######################  ##########################################################################
#                         ^
#                         |
####################### auxilliary functions ######################################################



####################### Wacker forward #############################################################
#                         |
#                         v 
####################################################################################################

"""
    Wackerforward(α_inc,  β_aut, Jθ, Jϕ[, Jχ])


Return regularly sampled S12 measurement data on sphere.
"""
function Wackerforward(
    α_inc,
    β_aut::Array{C,1},
    Jθ::Integer,
    Jϕ::Integer,
) where {C<:Complex}
    # First-order probe

    L, Nθ, Lθ = _inputdimensions(α_inc, β_aut, Jθ)
    u = _sum_product_αβ(α_inc, β_aut, L; firstorder = true)
    v_ =
        _zeropaddingθ!(_expandθmodes(u, L; firstorder = true), L, Jθ, Lθ; firstorder = true)
    v = _zeropaddingϕ!(fft!(v_, 1)[1:Nθ, :, :], L, Jϕ, Nθ; firstorder = true)
    w = fft!(v, 2)
    return _χmodes_to_S12_Firstorder(w)
end

function Wackerforward(
    α_inc::Array{C,1},
    β_aut::Array{C,1},
    Jθ::Integer,
    Jϕ::Integer,
    Jχ::Integer,
) where {C<:Complex}
    # Higher-order probe

    L, Nθ, Lθ = _inputdimensions(α_inc, β_aut, Jθ)
    u = _sum_product_αβ(α_inc, β_aut, L; firstorder = false)
    v_ = _zeropaddingθ!(
        _expandθmodes(u, L; firstorder = false),
        L,
        Jθ,
        Lθ;
        firstorder = false,
    )
    v = _zeropaddingϕ!(fft!(v_, 1)[1:Nθ, :, :], L, Jϕ, Nθ; firstorder = false)
    w = _zeropaddingχ!(fft!(v, 2), L, Jχ, Jϕ, Nθ)
    return fft(w, 3)
end
function Wackerforward(
    αin::Array{C,1},
    β_aut::Array{C,1},
    θvec::Array{F,1},
    Jϕ::Integer,
) where {C<:Complex} where {F<:Number}
    # First-order probe, irregular sampling in θ

    nthet = length(θvec)
    s = 0
    if isodd(length(αin))
        s = 1
    else
        s = 2
    end
    # L is the maximum mode order
    L = Int(floor(sqrt((length(αin) - s) / 2 + 1)))

    u = _sum_product_αβ(αin, β_aut, L; firstorder = true)

    v_ = zeros(C, nthet, 2 * L + 1, 2)

    cosmθ = [2 * cos.(mθ * θvec) for mθ = 1:L]
    sinmθ = [complex(0, 2) * sin.(mθ * θvec) for mθ = 1:L]
    for ℓ = 1:L

        Δℓ = _Δℓ_mμ(ℓ)

        indices = ifftshift(1:(2ℓ+1))

        for m = (-ℓ):ℓ
            mind = mod(m, 2 * L + 1) + 1 #m+L+1
            Δmind = indices[mod(m, 2 * ℓ + 1)+1] #m+ℓ+1
            k = ℓ * (ℓ + 1) + m
            for μ = -1:2:1
                Δμind = -1
                μind = mod(3 + μ, 3) # μind=1 => μ=1;   μind=2 => μ=-1
                if μ == 1
                    Δμind = indices[2]
                elseif μ == -1
                    Δμind = indices[2*ℓ+1]
                end
                mθind = indices[1]
                ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]
                imfac = ((1.0im)^(μ - m) * u[k, μind])

                vview = view(v_, :, mind, μind)
                vview .+= imfac .* ΔΔ

                trig = isodd(μ + m) ? sinmθ : cosmθ
                for mθ = 1:ℓ
                    mθind = indices[mod(mθ, 2 * ℓ + 1)+1]
                    ΔΔ = Δℓ[mθind, Δμind] * Δℓ[mθind, Δmind]

                    vview .+= imfac .* ΔΔ .* trig[mθ]
                end


            end
        end
    end
    v = _zeropaddingϕ!(v_, L, Jϕ, nthet)

    w = fft!(v, 2)
    return _χmodes_to_S12_Firstorder(w)
end

function Wackerforward(
    αin::Array{C,1},
    β_aut::Array{C,1},
    θϕvec::Array{Tuple{N,F},1},
) where {C<:Complex,N<:Number,F<:Real}
    # First-order probe, irregular sampling in θ and ϕ
    w = zeros(C, length(θϕvec), 2)

    L, _, __ = _inputdimensions(αin, β_aut, 0)
    u = _sum_product_αβ(αin, β_aut, L; firstorder = true)
    v_ = _expandθmodes(u, L; firstorder = true)
    mvec = (fftfreq(2 * L + 1, 2 * L + 1))
    Threads.@threads for k in eachindex(θϕvec)
        for mθ = 1:(2L+1), m = 1:(2*L+1)
            fac = cis(-mvec[mθ] * θϕvec[k][1] - mvec[m] * θϕvec[k][2])
            w[k, 1] += v_[mθ, m, 1] * fac
            w[k, 2] += v_[mθ, m, 2] * fac
        end
    end

    S12 = Array{C}(undef, size(w))
    S12[:, 1] .= w[:, 1] .+ w[:, 2]
    S12[:, 2] .= 1im .* (w[:, 1] .- w[:, 2])
    return S12
end

function Wackerforward(
    αin::Array{C,1},
    β_aut::Array{C,1},
    θϕvec::Array{Tuple{N,F},1},
    k0Δz::Real,
) where {C<:Complex,N<:Number,F<:Real}
    # First-order probe, irregular sampling in θ and ϕ, spectrum is evaluated at some distance Δz with wavenumber k
    w = zeros(C, length(θϕvec), 2)

    L, _, __ = _inputdimensions(αin, β_aut, 0)
    u = _sum_product_αβ(αin, β_aut, L; firstorder = true)
    v_ = _expandθmodes(u, L; firstorder = true)
    mvec = (fftfreq(2 * L + 1, 2 * L + 1))
    Threads.@threads for k in eachindex(θϕvec)
        k0Δzcosθ = k0Δz * cos(θϕvec[k][1])
        for mθ = 1:(2L+1), m = 1:(2*L+1)
            fac = cis(-mvec[mθ] * θϕvec[k][1] - mvec[m] * θϕvec[k][2] - k0Δzcosθ)
            w[k, 1] += v_[mθ, m, 1] * fac
            w[k, 2] += v_[mθ, m, 2] * fac
        end
    end

    S12 = Array{C}(undef, size(w))
    S12[:, 1] .= view(w, :, 1) .+ view(w, :, 2)
    S12[:, 2] .= 1im .* (view(w, :, 1) - view(w, :, 2))
    return S12
end



####################### Wacker forward  adjoint #############################################################
#                         |
#                         v 
####################################################################################################

"""
    Wackerforward_ad(α_inc, S21, Jθ, Jϕ, Jχ)
  
  
Return adjoint operator to `Wackerforward(α_inc,  β_aut, Jθ, Jϕ, Jχ)`.
"""
function Wackerforward_ad(
    α_inc::Array{C,1},
    S12::Array{K,3},
    Jθ::Integer,
    Jϕ::Integer,
    Jχ::Integer,
) where {C<:Complex} where {K<:Complex}
    L, Nθ, Lθ = _inputdimensions(α_inc, α_inc, Jθ)

    W = bfft(S12, 3)
    w = _zeropaddingχ_ad(W, L, Jχ, Jϕ, Nθ)
    V = bfft(w, 2)
    v = _zeropaddingϕ_ad(V, L, Jϕ, Nθ; firstorder = false)
    V_ = bfft(v, 1)
    v_ = _zeropaddingθ_ad(V_, L, Jθ, Lθ; firstorder = false)
    u = _expandθmodes_ad(v_, L; firstorder = false)
    β_aut = _sum_product_αβ_ad(u, α_inc, L; firstorder = false)
    return β_aut
end




####################### Wacker inverse #############################################################
#                         |
#                         v 
####################################################################################################

"""
    Wacker(S12, αin, [θvec, weightvec])


Return (receiving) spherical mode coefficients from measured S12 data with regular grid points on sphere.
S12 is expected to have dimensions (L+1) × (2L+1) × 2 if no θvec or Jθ, JΦ, Jχ is given
"""
function Wacker(S12::Array{C,3}, αin::Array{T,1}) where {C<:Complex} where {T<:Complex}
    ## | Determine size of input data
    ## v   
    a, b, c = size(S12)
    L = a - 1
    b == 2 * L + 1 ? () : throw(DimensionMismatch("Second dimension of S12 must be 2L+1."))
    c == 2 ? () : throw(DimensionMismatch("Third dimension of S12 must be 2."))
    ## ^
    ## | Determine size of input data

    v = _χandφ_integral(S12, L)
    v .= ifft!(v, 1)
    u = _θintegral(v, L)
    return _receivecoeffs_2by2matrix(u, αin, L)
end

function Wacker(
    S12::Array{C,3},
    αin::Array{C,1},
    θvec::Array{F,1},
    weightvec::Array{K,1},
) where {C<:Complex} where {F<:Real} where {K<:Number}
    ## | Determine size of input data
    ## v   
    a, b, c = size(S12)
    L = a - 1
    b == 2 * L + 2 ? () :
    throw(DimensionMismatch("Second dimension of S12 must be >=2L+2."))
    c == 2 ? () : throw(DimensionMismatch("Third dimension of S12 must be 2."))
    ## ^
    ## | Determine size of input data

    w = Array{C}(undef, a, b, 2)
    w[:, :, 1] = S12[:, :, 1] - 1im * S12[:, :, 2] # w_1(ϑ, ϕ)
    w[:, :, 2] = S12[:, :, 1] + 1im * S12[:, :, 2] # w_-1 (ϑ, ϕ)
    w .*= 0.5

    v_ = ifft!(w, 2)
    v = [v_[:, 1:(L+1), :] v_[:, (end-L+1):end, :]] # remove redundant sample of ϕ-harmonics

    u = zeros(C, L * (L + 2), 2)
    nthet = length(θvec)
    cosmθ = [2 * cos.(mθ * θvec) for mθ = 1:L]
    sinmθ = [complex(0, 2) * sin.(mθ * θvec) for mθ = 1:L]
    dvec = Array{ComplexF64}(undef, nthet)
    for ℓ = 1:L
        Δℓ = _Δℓ_mμ(ℓ)
        indices = ifftshift(1:(2ℓ+1))

        for m = (-ℓ):ℓ
            m_ind = mod(m, 2 * L + 1) + 1
            Δm_ind = indices[mod(m, 2 * ℓ + 1)+1]
            k = ℓ * (ℓ + 1) + m
            for μ = -1:2:1
                Δμ_ind = -1
                μ_ind = mod(3 + μ, 3)
                if μ == 1
                    Δμ_ind = indices[2]
                elseif μ == -1
                    Δμ_ind = indices[2*ℓ+1]
                end
                mθind = indices[1]
                ΔΔ = Δℓ[mθind, Δμ_ind] * Δℓ[mθind, Δm_ind]
                dvec .= ΔΔ
                trig = isodd(μ + m) ? sinmθ : cosmθ
                for mθ = 1:ℓ
                    mθind = indices[mod(mθ, 2 * ℓ + 1)+1]
                    ΔΔ = Δℓ[mθind, Δμ_ind] * Δℓ[mθind, Δm_ind]
                    dvec .+= ΔΔ .* trig[mθ]
                end
                # dvec .*= (1.0im)^(μ-m)

                # u[k, μ_ind] += (2ℓ + 1) / 2 .* sum(dvec .* weightvec .* v[:, m_ind, μ_ind])
                dvec .*= (1.0im)^(μ - m) .* weightvec
                u[k, μ_ind] += (2ℓ + 1) ./ 2 .* (transpose(dvec) * view(v, :, m_ind, μ_ind))
            end

        end

    end
    return _receivecoeffs_2by2matrix(u, αin, L)
end
