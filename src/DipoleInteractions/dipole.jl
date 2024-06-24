"""
Abstract type container type for HertzDipole and FitzgeraldDipole
"""
abstract type AbstractDipole <: AntennaFieldRepresentation end

"""
    HertzDipole{T<:Number} <: AbstractDipole
Store information for an electric dipole

# Fields 
- `pos::SArray{Tuple{3},T,1,3}` : position (can be real or complex)
- `dir::SArray{Tuple{3},T,1,3}` : direction (must be complex to enable circular polarization)
- `mag::ComplexF64` : magnitude (for scalar excitation)
"""
mutable struct HertzDipole{T,C} <: AbstractDipole where {T<:Number,C<:Complex}
    pos::SArray{Tuple{3},T,1,3}
    dir::SArray{Tuple{3},C,1,3}
    mag::C
end
# mutable struct HertzDipole{C} <: AbstractDipole where {C<:Complex} 
#     pos::SArray{Tuple{3},C,1,3}
#     dir::SArray{Tuple{3},C,1,3}
#     mag::C
# end
function HertzDipole(pos, dir, mag)
    return HertzDipole(SVector{3}(pos), SVector{3}(complex.(dir)), complex(mag))
end

function elementtype(dipole::HertzDipole{T,C}) where {T,C}
    return C
end

function elementtype(dipoles::Array{<:AbstractDipole})
    T = elementtype(dipoles[1])
    for dipole in dipoles
        T = promote_type(T, elementtype(dipole))
    end
    return T
end

"""
    FitzgeraldDipole{T<:Number} <: AbstractDipole
Store information for a magnetic dipole
 # Fields 
 - `pos::SArray{Tuple{3},T,1,3}` : position (can be real or complex)
 - `dir::SArray{Tuple{3},T,1,3}` : direction (must be complex to enable circular polarization)
 - `mag::ComplexF64` : magnitude (for scalar excitation)
"""
mutable struct FitzgeraldDipole{T,C} <: AbstractDipole where {T<:Number,C<:Complex}
    pos::SArray{Tuple{3},T,1,3}
    dir::SArray{Tuple{3},C,1,3}
    mag::C
end
# mutable struct FitzgeraldDipole{C} <: AbstractDipole where {C<:Complex} 
#     pos::SArray{Tuple{3},C,1,3}
#     dir::SArray{Tuple{3},C,1,3}
#     mag::C
# end
function FitzgeraldDipole(pos, dir, mag)
    return FitzgeraldDipole(SVector{3}(pos), SVector{3}(convert.(ComplexF64, dir)), convert(ComplexF64, mag))
end
function elementtype(dipole::FitzgeraldDipole{T,C}) where {T,C}
    return C
end


#### Type Conversions
import Base.complex
function complex(dipole::HertzDipole{<:Number})
    #type conversion for dipoles with real position to complex position
    pos = SVector(complex(dipole.pos))
    return HertzDipole(pos, dipole.dir, dipole.mag)
end
function complex(diparray::Array{HertzDipole{T},1}) where {T<:Number}
    #type conversion for dipoles with real position to complex position
    newarray = Array{HertzDipole{complex(T)},1}(undef, length(diparray))
    for k in eachindex(diparray)
        newarray[k] = complex(diparray[k])
    end
    return newarray
end

function converttype(T::Type{FitzgeraldDipole{N}}, dipole::HertzDipole{<:Number}) where {N<:Number}
    return FitzgeraldDipole(convert.(N, dipole.pos), dipole.dir, dipole.mag)
end
function converttype(T::HertzDipole{N}, dipole::FitzgeraldDipole{<:Number}) where {N<:Number}
    return HertzDipole(convert.(N, dipole.pos), dipole.dir, dipole.mag)
end


"""
    g₀(R, k₀) -> Number
Scalar Greens function exp(-j*k₀*R)/(4*π*R)
"""
function g₀(R::Array{<:Number,3}, k₀::Float64)
    r = cdist(R)
    return cis(-k₀ * r) / (4 * π * r) # exp(-j*k₀*r)/(4*π*r) 
end
function g₀(R::SArray{Tuple{3},C,1,3}, k₀::Float64) where {C<:Number}
    r = cdist(R)
    return cis(-k₀ * r) / (4 * π * r) # exp(-j*k₀*r)/(4*π*r) 
end
function g₀(R::Number, k₀::Float64)
    return cis(-k₀ * R) / (4 * π * R) # exp(-j*k₀*r)/(4*π*r) 
end


"""
    Greensdivdiv(R::SArray{Tuple{3},C,1,3}, dipoledir, k₀::Float64) -> SArray{Tuple{3},C,1,3}
    Greensdivdiv(R::Array{<:Number,3}, dipoledir, k₀::Float64) -> Array{<:Number,3}

Return  dyadic operator (I+1/k² ∇∇) g₀(r,r') 
"""
function Greensdivdiv(R, dipoledir, k₀::Float64)

    d = cdist(R)
    kd = (k₀ * d)
    kd² = kd^2
    ℓeᵣeᵣ = R * (udot(R, dipoledir)) / (d^2)

    return (g₀(d, k₀) * ((3.0 / kd² + 3.0im / kd - 1) * ℓeᵣeᵣ + (-1.0 / kd² - 1.0im / kd + 1) * dipoledir)[:])
end


"""
  Greensrot(R::SArray{Tuple{3},C,1,3}, dipoledir, k₀::Float64) -> SArray{Tuple{3},C,1,3}
  Greensrot(R::Array{<:Number,3}, dipoledir, k₀::Float64) -> Array{<:Number,3}

Return dyadic operator ∇g₀(r,r') × I
"""
function Greensrot(R, dipoledir, k₀::Float64)

    d = cdist(R)
    eᵣ = real(R) / norm(real(R))
    return cross(eᵣ, dipoledir) * g₀(d, k₀) * (-1im * k₀ - 1 / d)
end


function transmission(
    sourcedipole::HertzDipole{T,C1}, probedipole::HertzDipole{F,C2}, k₀::Float64
) where {T<:Number,F<:Number,C1<:Complex,C2<:Complex}
    C = promote_type(C1, C2)
    # HertzDipole receives Efield
    R = probedipole.pos - sourcedipole.pos
    vec = Greensdivdiv(R, sourcedipole.dir, k₀)
    b = udot(vec, probedipole.dir) #)vec[1] * probedipole.dir[1] + vec[2] * probedipole.dir[2] + vec[3] * probedipole.dir[3]
    return convert(C, -0.5 * 1im * k₀ * Z₀ * b * sourcedipole.mag * probedipole.mag)
end
function transmission(
    sourcedipole::FitzgeraldDipole{T,C1}, probedipole::FitzgeraldDipole{F,C2}, k₀::Float64
) where {T<:Number,F<:Number,C1<:Complex,C2<:Complex}
    C = promote_type(C1, C2)
    # FitzgeraldDipole receives negative Hfield
    R = probedipole.pos - sourcedipole.pos
    vec = Greensdivdiv(R, sourcedipole.dir, k₀)
    b = udot(vec, probedipole.dir) #)vec[1] * probedipole.dir[1] + vec[2] * probedipole.dir[2] + vec[3] * probedipole.dir[3]
    return convert(C, 0.5 * 1im * k₀ / Z₀ * b * sourcedipole.mag * probedipole.mag) #missing - sign because Fitzgerald dipole receives negative Hfield
end
function transmission(
    sourcedipole::FitzgeraldDipole{T,C1}, probedipole::HertzDipole{F,C2}, k₀::Float64
) where {T<:Number,F<:Number,C1<:Complex,C2<:Complex}
    C = promote_type(C1, C2)
    # HertzDipole receives Efield
    R = probedipole.pos - sourcedipole.pos
    vec = Greensrot(R, sourcedipole.dir, k₀)
    b = udot(vec, probedipole.dir) #)vec[1] * probedipole.dir[1] + vec[2] * probedipole.dir[2] + vec[3] * probedipole.dir[3]
    return convert(C, -0.5 * b * sourcedipole.mag * probedipole.mag)
end
function transmission(
    sourcedipole::HertzDipole{T,C1}, probedipole::FitzgeraldDipole{F,C2}, k₀::Float64
) where {T<:Number,F<:Number,C1<:Complex,C2<:Complex}

    # FitzgeraldDipole receives negative Hfield
    return transmission(probedipole, sourcedipole, k₀)
end
function transmission(field, dipoles::Array{<:AbstractDipole}, k₀::Number)
    R = [dipole.pos for dipole in dipoles]
    fieldtuples = ehfield(field, R, k₀)
    E = [fieldtuple[1] for fieldtuple in fieldtuples]
    H = [fieldtuple[2] for fieldtuple in fieldtuples]

    T = promote_type(elementtype(field), elementtype(dipoles))

    b = zero(T)
    for k in eachindex(dipoles)
        if isa(dipoles[k], HertzDipole{<:Number})
            b += 0.5 * sum(E[k] .* dipoles[k].dir) * dipoles[k].mag
        elseif isa(dipoles[k], FitzgeraldDipole{<:Number})
            b += -0.5 * udot(H[k], dipoles[k].dir) * dipoles[k].mag
        end
    end
    return b
end



function dipolefarfield(
    sourcedipole::HertzDipole{T,C}, eᵣ::SArray{Tuple{3},N,1,3}, eθ::SArray{Tuple{3},N,1,3}, eϕ::SArray{Tuple{3},N,1,3}, k₀::Float64
) where {T<:Number,N<:Number,C<:Complex}

    E_FF = C(0.0, -k₀) * Z₀ / (4π) * sourcedipole.mag * cis(k₀ * udot(eᵣ, sourcedipole.pos))
    Eθ = C(udot(eθ, sourcedipole.dir) * E_FF)
    Ephi = C(udot(eϕ, sourcedipole.dir) * E_FF)
    return Eθ, Ephi
end
function dipolefarfield(
    sourcedipole::FitzgeraldDipole{T,C},
    eᵣ::SArray{Tuple{3},N,1,3},
    eθ::SArray{Tuple{3},N,1,3},
    eϕ::SArray{Tuple{3},N,1,3},
    k₀::Float64,
) where {T<:Number,N<:Number,C<:Complex}

    #far-field of Fitzgerald dipole is scaled by 1/Z₀ and rotated
    E_FF = C(0.0, -k₀) / (4π) * sourcedipole.mag * cis(k₀ * udot(eᵣ, sourcedipole.pos))
    Eθ = C(udot(eϕ, sourcedipole.dir) * E_FF)
    Ephi = C(-udot(eθ, sourcedipole.dir) * E_FF)
    return Eθ, Ephi
end
function farfield(sourcedipole::AbstractDipole, θ::Number, ϕ::Number, k₀::Float64)
    st, ct = sincos(θ)
    sp, cp = sincos(convert(typeof(θ), ϕ))


    eᵣ = SVector(st * cp, st * sp, ct)
    eθ = SVector(ct * cp, ct * sp, -st)
    eϕ = SVector(-sp, cp, 0.0)

    return dipolefarfield(sourcedipole, eᵣ, eθ, eϕ, k₀)
end

function farfield(sourcedipoles::Array{D,1}, θs::Array{<:Number,1}, ϕs::Array{<:Number,1}, k₀::Float64) where{D<:Union{HertzDipole{T,C}, FitzgeraldDipole{T,C}} where{T,C}}
    C = elementtype(sourcedipoles)
    # Eθ = zeros(C, length(θs), length(ϕs))
    # Eϕ = zeros(C, length(θs), length(ϕs))
    Eθ = Matrix{C}(undef, length(θs), length(ϕs))
    Eϕ = Matrix{C}(undef, length(θs), length(ϕs))

    sp = sin.(convert.(eltype(θs), ϕs))
    cp = cos.(convert.(eltype(θs), ϕs))
    st = sin.(θs)
    ct = cos.(θs)
    for kϕ in eachindex(ϕs)
        eϕ = SVector(-sp[kϕ], cp[kϕ], 0.0) 
        
            for kθ in eachindex(θs)
            

            eᵣ = SVector(st[kθ] * cp[kϕ], st[kθ] * sp[kϕ], ct[kθ])
            eθ = SVector(ct[kθ] * cp[kϕ], ct[kθ] * sp[kϕ], -st[kθ])
            reset= true
            for kd in eachindex(sourcedipoles)
                Etheta, Ephi = dipolefarfield(sourcedipoles[kd], eᵣ, eθ, eϕ, k₀)
                if reset 
                    Eθ[kθ, kϕ] = Etheta
                    Eϕ[kθ, kϕ] = Ephi
                else
                    Eθ[kθ, kϕ] += Etheta
                    Eϕ[kθ, kϕ] += Ephi
                    reset =false
                end
            end
        end
    end
    return Eθ, Eϕ
end


function efield(dipole::AbstractDipole, R::AbstractVector{<:Number}, k₀::Number)
    C = elementtype(dipole)
    xdipole = HertzDipole(R, [1; 0; 0], convert(C, (1.0)))
    ydipole = HertzDipole(R, [0; 1; 0], convert(C, (1.0)))
    zdipole = HertzDipole(R, [0; 0; 1], convert(C, (1.0)))
    return SVector{3}(2 * [transmission(dipole, xdipole, k₀); transmission(dipole, ydipole, k₀); transmission(dipole, zdipole, k₀)])
end
function efield(dipoles::Array{<:AbstractDipole}, R::AbstractVector{<:Number}, k₀::Real)
    C = elementtype(dipoles)
    xdipole = HertzDipole(R, [1; 0; 0], convert(C, (1.0)))
    ydipole = HertzDipole(R, [0; 1; 0], convert(C, (1.0)))
    zdipole = HertzDipole(R, [0; 0; 1], convert(C, (1.0)))
    Ex = zero(C)
    Ey = zero(C)
    Ez = zero(C)
    for dipole in dipoles
        Ex += 2 * transmission(dipole, xdipole, k₀)
        Ey += 2 * transmission(dipole, ydipole, k₀)
        Ez += 2 * transmission(dipole, zdipole, k₀)
    end
    return SVector{3}([Ex, Ey, Ez])
end



function hfield(dipole::AbstractDipole, R::AbstractVector{<:Number}, k₀::Number)
    C = elementtype(dipole)
    xdipole = FitzgeraldDipole(R, [1; 0; 0], convert(C, (1.0)))
    ydipole = FitzgeraldDipole(R, [0; 1; 0], convert(C, (1.0)))
    zdipole = FitzgeraldDipole(R, [0; 0; 1], convert(C, (1.0)))
    return SVector{3}(-2 * [transmission(dipole, xdipole, k₀); transmission(dipole, ydipole, k₀); transmission(dipole, zdipole, k₀)])
end
function hfield(dipoles::Array{<:AbstractDipole}, R::AbstractVector{<:Number}, k₀::Real)
    C = elementtype(dipoles)
    xdipole = FitzgeraldDipole(R, [1; 0; 0], convert(C, (1.0)))
    ydipole = FitzgeraldDipole(R, [0; 1; 0], convert(C, (1.0)))
    zdipole = FitzgeraldDipole(R, [0; 0; 1], convert(C, (1.0)))
    Hx = zero(C)
    Hy = zero(C)
    Hz = zero(C)
    for dipole in dipoles
        Hx += -2 * transmission(dipole, xdipole, k₀)
        Hy += -2 * transmission(dipole, ydipole, k₀)
        Hz += -2 * transmission(dipole, zdipole, k₀)
    end
    return SVector{3}([Hx, Hy, Hz])
end


function ehfield(dipole::AbstractDipole, R, k₀::Number)
    E = Efield(dipole, R, k₀)
    H = Hfield(dipole, R, k₀)
    return E, H
end
function ehfield(dipoles::Array{<:AbstractDipole,1}, R::AbstractVector{<:Number}, k₀::Real)
    E = sum([efield(dipole, R, k₀) for dipole in dipoles])
    H = sum([hfield(dipole, R, k₀) for dipole in dipoles])
    return E, H
end
function rotate(dipole::AbstractDipole, χ::Number, θ::Number, ϕ::Number)
    newdipole = deepcopy(dipole)
    R = rot_mat_zyz(χ, θ, ϕ)
    newdipole.pos = R * dipole.pos
    newdipole.dir = R * dipole.dir
    return newdipole
end
function rotate(dipole::AbstractDipole; χ=0.0, θ=0.0, ϕ=0.0)
    return rotate(dipole, χ, θ, ϕ)
end




