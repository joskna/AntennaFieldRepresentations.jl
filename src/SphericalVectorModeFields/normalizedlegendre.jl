"""
    sphPlm_deriv_array(Lmax,m,x)

Return normalized Associated Legendre polynomials and their derivaives up to Lmax
sphPlm= 1/√(2π)  √((2ℓ+1)(ℓ-m)!/(2(ℓ+m)!)  Plm

From [T. Limpanuparb, J. Milthorpe: "Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications", arXiv 2014](https://arxiv.org/abs/1410.1748)
"""
function sphPlm_deriv_array(Lmax::I, m::I, x::T) where {I<:Integer} where {T<:Number}
    ysquared = (1 - x^2)
    y = (√ysquared)
    msquared = m^2
    # sinx=sin(x)
    # cosx=cos(x)

    # use one-term recurrence relation for Pℓℓ in direction of increasing ℓ 
    # Pℓℓ=0.70710678118654757273731092936941422522068023681640625 # √0.5
    Pℓℓ = 0.282094791773878139640174822488916106522083282470703125 # 0.5*√(1/π)
    # Pℓℓ=0.56418958354775627928034964497783221304416656494140625 #√(1/π)
    for ℓ = 1:m
        Pℓℓ = -√(1 + 1 / (2 * ℓ)) * y * Pℓℓ
    end

    # use two-term recurrence relation for Pℓm in direction of increasing ℓ 
    sphPlm = zeros(T, Lmax - m + 1)
    sphPlmderiv = zeros(T, Lmax - m + 1)
    sphPlm[1] = Pℓℓ
    sphPlmderiv[1] = -m * x / ysquared * sphPlm[1]
    if Lmax > m
        ℓ = m + 1
        ℓsquared = ℓ^2
        sphPlm[2] = √(2m + 3) * x * sphPlm[1]
        sphPlmderiv[2] =
            (
                -ℓ * x * sphPlm[2] +
                √((2 * ℓ + 1) * (ℓsquared - msquared) / (2 * ℓ - 1)) * sphPlm[1]
            ) / ysquared


        for k = 3:(Lmax-m+1)
            ℓ = (m + k - 1)
            ℓsquared = ℓ^2
            ℓm1squared = (ℓ - 1)^2
            aℓm = √((4 * ℓsquared - 1) / (ℓsquared - msquared))
            bℓm = -√((ℓm1squared - msquared) / (4 * ℓm1squared - 1))
            sphPlm[k] = aℓm * (x * sphPlm[k-1] + bℓm * sphPlm[k-2])
            sphPlmderiv[k] =
                (
                    -ℓ * x * sphPlm[k] +
                    √((2 * ℓ + 1) * (ℓsquared - msquared) / (2 * ℓ - 1)) * sphPlm[k-1]
                ) / ysquared
        end
    end

    return sphPlm, sphPlmderiv
end

"""
    Plm_deriv_array(Lmax,m,x)

Return unnormalized Associated Legendre polynomials and their derivaives up to Lmax
"""
function Plm_deriv_array(Lmax::I, m::I, x::T) where {I<:Integer} where {T<:Number}
    ysquared = (1 - x^2)
    y = √ysquared

    # sinx=sin(x)
    # cosx=cos(x)

    # use one-term recurrence relation for Pℓℓ in direction of increasing ℓ 
    # Pℓℓ=0.70710678118654757273731092936941422522068023681640625 # √0.5
    # Pℓℓ=0.282094791773878139640174822488916106522083282470703125 # 0.5*√(1/π)
    Pℓℓ = one(T)
    for ℓ = 1:m
        Pℓℓ = -(2 * ℓ - 1) * y * Pℓℓ
    end

    # use two-term recurrence relation for Pℓm in direction of increasing ℓ 
    Plm = zeros(T, Lmax - m + 1)
    Plmderiv = zeros(T, Lmax - m + 1)
    Plm[1] = Pℓℓ
    Plmderiv[1] = -m * x / ysquared * Plm[1]
    if Lmax > m
        ℓ = m
        Plm[2] = x * (2 * ℓ + 1) * Plm[1]
        Plmderiv[2] = (ℓ * x * Plm[2] - (ℓ + m) * Plm[2]) / ysquared

        for k = 2:(Lmax-m)
            ℓ = (m + k - 1)
            Plm[k+1] = ((2ℓ + 1)x * Plm[k] - (ℓ + m) * Plm[k-1]) / (ℓ - m + 1)
            Plmderiv[k+1] = ((ℓ) * x * Plm[k+1] - (ℓ + m) * Plm[k]) / ysquared
        end
    end

    return Plm, Plmderiv
end

"""
    legendre_deps_array(m , Lmax, ϑ)

Return the theta dependencies occuring in the spherical wave expansion
"""
function legendre_deps_array(m::Integer, Lmax::Integer, ϑ::T) where {T<:Number}
    cost = cos(ϑ)
    sint = sin(ϑ)
    am = abs(m)

    Pℓm, Pℓmderiv = sphPlm_deriv_array(Lmax, am, cost)
    Pℓmderiv *= sint

    P3 = Pℓm
    P1 = m * P3 / sint
    P2 = -Pℓmderiv


    return P1, P2, P3
end

# """
#     Pl_array(Lmax,x)

# Return Legendre polynomials up to Lmax
# """
# function Pl_array(Lmax::I, x::T) where {I<:Integer} where {T<:Number}

#     Pℓ = zeros(T, Lmax + 1)
#     Pℓ[1] = one(T)
#     if Lmax > 0
#         Pℓ[2] = x
#     end

#     # use two-term recurrence relation for Pℓ in direction of increasing ℓ
#     for ℓ in 2:Lmax
#         Pℓ[ℓ + 1] = ((2 * ℓ - 1) * x * Pℓ[ℓ] - (ℓ - 1) * Pℓ[ℓ - 1]) / ℓ
#     end

#     return Pℓ
# end
