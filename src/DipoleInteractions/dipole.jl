


"""
    DipoleArray{P, E, T, C} <: AntennaFieldRepresentation{P, C}

Array of multiple small dipoles.

Behaves like an `AbstractVector{C}` with extra context.

# Type Parameters
- `P <: PropagationType`
- `E <: ElmagType`
- `T <: Number`
- `C <: Number`
"""
struct DipoleArray{P<:PropagationType,E<:ElmagType,T<:Number,C<:Number} <:
       AntennaFieldRepresentation{P,C}
    positions::Vector{SVector{3,T}}
    orientations::Vector{SVector{3,C}}
    dipolemoments::Vector{C}
    wavenumber::Number
    function DipoleArray{P,E,T,C}(
        positions::Vector{V1},
        orientations::Vector{V2},
        dipolemoments::Vector{C},
        wavenumber::Number,
    ) where {P<:PropagationType,E<:ElmagType,T<:Number,C<:Number,V1,V2}

        (length(orientations) != length(positions)) && throw(
            DimensionMismatch(
                "Input vectors `positions` and `orientations` must have the same lengths.",
            ),
        )
        (length(dipolemoments) != length(positions)) && throw(
            DimensionMismatch(
                "Input vectors `positions` and `dipolemoments` must have the same lengths.",
            ),
        )

        return new{P,E,T,C}(
            SVector{3,T}.(positions),
            SVector{3,C}.(orientations),
            dipolemoments,
            wavenumber,
        )
    end
end
function DipoleArray{P,E}(
    positions::Vector{V1},
    orientations::Vector{V2},
    dipolemoments::Vector{C},
    wavenumber,
) where {P<:PropagationType,E<:ElmagType,C,V1<:AbstractVector,V2<:AbstractVector{C}}
    return DipoleArray{P,E,eltype(V1),C}(positions, orientations, dipolemoments, wavenumber)
end

"""
    HertzArray{T, C}

Alias for `DipoleArray{Radiated,Electric,T,C}`

A `HertzArray` is used to represent an array of multiple small radiating electric dipoles.    
Behaves like an `AbstractVector{C}` with extra context.

# Type Parameters
- `T <: Number`
- `C <: Number`
"""
HertzArray{T,C} = DipoleArray{Radiated,Electric,T,C}
function HertzArray(
    positions::Vector{V1},
    orientations::Vector{V2},
    dipolemoments::Vector{C},
    wavenumber,
) where {C,V1<:AbstractVector,V2<:AbstractVector{C}}
    return HertzArray{eltype(V1),C}(positions, orientations, dipolemoments, wavenumber)
end

"""
    FitzgeraldArray{T, C}

Alias for `DipoleArray{Radiated,Magnetic,T,C}`

A `FitzgeraldArray` is used to represent an array of multiple small radiating magnetic dipoles.    
Behaves like an `AbstractVector{C}` with extra context.

# Type Parameters
- `T <: Number`
- `C <: Number`
"""
FitzgeraldArray{T,C} = DipoleArray{Radiated,Magnetic,T,C}
function asvector(dips::DipoleArray)
    return dips.dipolemoments
end
function FitzgeraldArray(
    positions::Vector{V1},
    orientations::Vector{V2},
    dipolemoments::Vector{C},
    wavenumber,
) where {C,T,V1<:AbstractVector{T},V2<:AbstractVector{C}}
    return FitzgeraldArray{T,C}(positions, orientations, dipolemoments, wavenumber)
end
function Base.similar(dips::DipoleArray{P,E,T,C}) where {P,E,T,C}
    return DipoleArray{P,E,T,C}(
        deepcopy(dips.positions),
        deepcopy(dips.orientations),
        similar(dips.dipolemoments),
        dips.wavenumber,
    )
end
Base.size(dips::DipoleArray) = size(dips.dipolemoments)
Base.getindex(dips::DipoleArray, i) = getindex(dips.dipolemoments, i)
Base.setindex!(dips::DipoleArray, v, i) = (dips.dipolemoments, v, i)
function setwavenumber!(dips::DipoleArray{P,E,C,T}, val::T) where {P,E,C,T}
    dips = DipoleArray{P,E,C,T}(dips.positions, dips.orientations, dips.dipolemoments, val)
    return dips
end

"""
    g₀(R::Number, k₀::Number) -> Number

Scalar Greens function exp(-j*k₀*R)/(4*π*R)
"""
function g₀(R::Number, k₀::Number)
    return g₀(R, k₀, Radiated)
end
function g₀(R::Number, k₀::Number, ::Radiated)
    return cis(-k₀ * R) / (4 * π * R) # exp(-j*k₀*r)/(4*π*r) 
end
function g₀(R::Number, k₀::Number, ::Absorbed)
    return cis(k₀ * R) / (4 * π * R) # exp(-j*k₀*r)/(4*π*r) 
end
function g₀(R::Number, k₀::Number, ::Incident)
    return sinc(k₀ * R / π) / (4 * π) # sin(k₀*r)/(4*π*r) 
end


"""
    Greensdivdiv(R, dipoledir, k₀)

Return  dyadic operator (I+1/k² ∇∇) g₀(r,r') 
"""
function Greensdivdiv(R, dipoledir, k₀)
    return Greensdivdiv(R, dipoledir, k₀, Radiated())
end
function Greensdivdiv(R, dipoledir, k₀, ::Radiated)

    d = cdist(R)
    kd = (k₀ * d)
    kd² = kd^2
    ℓeᵣeᵣ = R * (udot(R, dipoledir)) / (d^2)

    return (
        g₀(d, k₀, Radiated()) *
        ((3.0/kd²+3.0im/kd-1)*ℓeᵣeᵣ+(-1.0/kd²-1.0im/kd+1)*dipoledir)[:]
    )
end
function Greensdivdiv(R, dipoledir, k₀, ::Absorbed)

    d = cdist(R)
    kd = (k₀ * d)
    kd² = kd^2
    ℓeᵣeᵣ = R * (udot(R, dipoledir)) / (d^2)

    return (
        g₀(d, k₀, Absorbed()) *
        ((3.0/kd²-3.0im/kd-1)*ℓeᵣeᵣ+(-1.0/kd²+1.0im/kd+1)*dipoledir)[:]
    )
end
function Greensdivdiv(R, dipoledir, k₀, ::Incident)
    return (
        Greensdivdiv(R, dipoledir, k₀, Absorbed()) -
        Greensdivdiv(R, dipoledir, k₀, Radiated())
    ) ./ 2
end


"""
  Greensrot(R, dipoledir, k₀)

Return dyadic operator ∇g₀(r,r') × I
"""
function Greensrot(R, dipoledir, k₀)
    return Greensrot(R, dipoledir, k₀, Radiated())
end
function Greensrot(R, dipoledir, k₀, ::Radiated)

    d = cdist(R)
    eᵣ = real(R) / norm(real(R))
    return cross(eᵣ, dipoledir) * g₀(d, k₀, Radiated()) * (-1im * k₀ - 1 / d)
end
function Greensrot(R, dipoledir, k₀, ::Absorbed)

    d = cdist(R)
    eᵣ = real(R) / norm(real(R))
    return cross(eᵣ, dipoledir) * g₀(d, k₀, Absorbed) * (1im * k₀ - 1 / d)
end
function Greensrot(R, dipoledir, k₀, ::Incident)
    return (
        Greensrot(R, dipoledir, k₀, Absorbed()) - Greensrot(R, dipoledir, k₀, Radiated())
    ) / 2
end

function farfield(
    dipoles::DipoleArray{Radiated,E,T,C},
    θϕ::Tuple{R,R},
) where {C,E<:ElmagType,T,R<:Real}
    θ, ϕ = θϕ
    st, ct = sincos(θ)
    sp, cp = sincos(convert(typeof(θ), ϕ))

    eᵣ = SVector(st * cp, st * sp, ct)
    eθ = SVector(ct * cp, ct * sp, -st)
    eϕ = SVector(-sp, cp, 0.0)
    Eθ = zero(C)
    Eϕ = zero(C)

    return _dipolefarfield!(dipoles, eᵣ, eθ, eϕ, Eθ, Eϕ, reset = false)
end

function _dipolefarfield!(
    dipoles::DipoleArray{Radiated,E,T,C},
    eᵣ::SArray{Tuple{3},T,1,3},
    eθ::SArray{Tuple{3},T,1,3},
    eϕ::SArray{Tuple{3},T,1,3},
    Eθ::C,
    Eϕ::C;
    reset = true,
) where {C,E<:ElmagType,T}

    reset && (Eθ = zero(C))
    reset && (Eϕ = zero(C))

    k₀ = dipoles.wavenumber


    E_FF =
        C(0.0, -k₀) * _dipolefarfieldscalingfactor(E()) / (4π) * dipoles.dipolemoments .*
        cis.(k₀ * udot.(Ref(eᵣ), dipoles.positions))

    for (i, dir) in enumerate(dipoles.orientations)
        Epolθ, Epolϕ = _dipoledarfieldpolarization(eθ, eϕ, dir, E())
        Eθ += E_FF[i] * Epolθ
        Eϕ += E_FF[i] * Epolϕ
    end

    return Eθ, Eϕ
end

function _dipoledarfieldpolarization(eθ, eϕ, dir::SVector{3,C}, ::Electric) where {C}
    Eθ = C(udot(eθ, dir))
    Eϕ = C(udot(eϕ, dir))
    return Eθ, Eϕ
end
function _dipolefarfieldpolarization(eθ, eϕ, dir::SVector{3,C}, ::Magnetic) where {C}
    Eθ = C(udot(eϕ, dir))
    Eϕ = C(-udot(eθ, dir))
    return Eθ, Eϕ
end
function _dipolefarfieldscalingfactor(::Electric)
    return Z₀
end
function _dipolefarfieldscalingfactor(::Magnetic)
    return 1
end

function GreensE(R, dipoledir, k₀, ::Electric)
    return GreensE(R, dipoledir, k₀, Electric(), Radiated())
end
function GreensE(R, dipoledir, k₀, ::Electric, ::Radiated)
    return -1im * k₀ * Z₀ * Greensdivdiv(R, dipoledir, k₀, Radiated())
end
function GreensE(R, dipoledir, k₀, ::Electric, ::Absorbed)
    return 1im * k₀ * Z₀ * Greensdivdiv(R, dipoledir, k₀, Absorbed())
end
function GreensE(R, dipoledir, k₀, ::Electric, ::Incident)
    return (
        GreensE(R, dipoledir, k₀, Electric(), Absorbed()) -
        GreensE(R, dipoledir, k₀, Electric(), Radiated())
    ) / 2
end
# TODO: remove redundant code with parametric P instead of Radiated, Absorbed etc.
function GreensE(R, dipoledir, k₀, ::Magnetic)
    return GreensE(R, dipoledir, k₀, Magnetic(), Radiated())
end
function GreensE(R, dipoledir, k₀, ::Magnetic, ::Radiated)
    return -Greensrot(R, dipoledir, k₀, Radiated())
end
function GreensE(R, dipoledir, k₀, ::Magnetic, ::Absorbed)
    return -Greensrot(R, dipoledir, k₀, Absorbed())
end
function GreensE(R, dipoledir, k₀, ::Magnetic, ::Incident)
    return -Greensrot(R, dipoledir, k₀, Incident())
end
function GreensH(R, dipoledir, k₀, ::Electric)
    return GreensH(R, dipoledir, k₀, Electric(), Radiated())
end
function GreensH(R, dipoledir, k₀, ::Electric, ::Radiated)
    return Greensrot(R, dipoledir, k₀, Radiated())
end
function GreensH(R, dipoledir, k₀, ::Electric, ::Absorbed)
    return Greensrot(R, dipoledir, k₀, Absorbed())
end
function GreensH(R, dipoledir, k₀, ::Electric, ::Incident)
    return Greensrot(R, dipoledir, k₀, Incident())
end
function GreensH(R, dipoledir, k₀, ::Magnetic)
    return GreensH(R, dipoledir, k₀, Magnetic(), Radiated())
end
function GreensH(R, dipoledir, k₀, ::Magnetic, ::Radiated)
    return -1im * k₀ / Z₀ * Greensdivdiv(R, dipoledir, k₀, Radiated())
end
function GreensH(R, dipoledir, k₀, ::Magnetic, ::Absorbed)
    return 1im * k₀ / Z₀ * Greensdivdiv(R, dipoledir, k₀, Absorbed())
end
function GreensH(R, dipoledir, k₀, ::Magnetic, ::Incident)
    return (
        GreensH(R, dipoledir, k₀, Magnetic(), Absorbed()) -
        GreensH(R, dipoledir, k₀, Magnetic(), Radiated())
    ) / 2
end

function efield!(
    storage,
    dipoles::DipoleArray{P,E,T,C},
    R;
    reset = true,
) where {P<:PropagationType,C,E<:ElmagType,T}

    k₀ = dipoles.wavenumber

    reset && fill!(storage, zero(C))

    for k in eachindex(dipoles.positions)
        Rvec = R - dipoles.positions[k]
        storage .+=
            GreensE(Rvec, dipoles.orientations[k], k₀, E(), P()) .* dipoles.dipolemoments[k]
    end
    return storage
end
function hfield!(
    storage,
    dipoles::DipoleArray{P,E,T,C},
    R;
    reset = true,
) where {P<:PropagationType,C,E<:ElmagType,T}

    k₀ = dipoles.wavenumber

    reset && fill!(storage, zero(C))

    for k in eachindex(dipoles.positions)
        Rvec = R - dipoles.positions[k]
        storage .+=
            GreensH(Rvec, dipoles.orientations[k], k₀, E(), P()) .* dipoles.dipolemoments[k]
    end
    return storage
end

function rotate!(
    rotated_dipoles::DipoleArray{P,E,T,C},
    dipoles::DipoleArray{P,E,T,C},
    χ,
    θ,
    ϕ,
) where {P<:PropagationType,C,E<:ElmagType,T}

    R = rot_mat_zyz(χ, θ, ϕ)
    for k in eachindex(dipoles.positions)
        rotated_dipoles.positions[k] = SVector{3}(R * dipoles.positions[k])
        rotated_dipoles.orientations[k] = SVector{3}(R * dipoles.orientations[k])
        rotated_dipoles.dipolemoments[k] = dipoles.dipolemoments[k]
    end
    return rotated_dipoles
end

function spatialshift!(
    shifted_dipoles::DipoleArray{P,E,T,C},
    dipoles::DipoleArray{P,E,T,C},
    R,
) where {P<:PropagationType,C,E<:ElmagType,T}
    for k in eachindex(shifted_dipoles.positions)
        shifted_dipoles.positions[k] .= dipoles.positions[k] .+ R
        shifted_dipoles.orientations[k] .= dipoles.orientations[k]
        shifted_dipoles.dipolemoments[k] = dipoles.dipolemoments
    end
    return shifted_dipoles
end

function _rmax(dipoles::DipoleArray)
    rmax = 0.0
    for pos in dipoles.positions
        if norm(pos) > rmax
            rmax = norm(pos)
        end
    end
    return rmax
end

function _rmin(dipoles::DipoleArray)
    rmin = Inf
    for pos in dipoles.positions
        if norm(pos) < rmin
            rmin = norm(pos)
        end
    end
    return rmin
end
function equivalentorder(dipoles::DipoleArray{P,E,T,C}; ϵ = 1e-7) where {P,E,T,C}
    k0 = getwavenumber(dipoles)
    rsph = (P() == Radiated()) ? (2 * _rmax(dipoles)) : (2 * _rmin(dipoles))
    L = _modeorder(rsph, k0; ϵ = ϵ)
    return L
end

function _modeorder(rmax::Real, wavenumber::Real; ϵ = 1e-7)
    L = maximum([
        3,
        Int(
            ceil(
                wavenumber * rmax +
                1.8 * (log10(1 / ϵ))^(2 / 3) * (wavenumber * rmax)^(1 / 3),
            ),
        ),
    ])
    return L
end
