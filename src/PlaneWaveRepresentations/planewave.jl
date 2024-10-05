
"""
    PlaneWaveExpansion{P <: PropagationType, S <: SphereSamplingStrategy, C} <: AntennaFieldRepresentation{P, C}

Collection of electromagnetic plane waves propagating into various directions.

Behaves like an `AbstractVector{C}` with extra context.
The sampling strategy stored by the type parameter `S` defines the propagation directions of the plane waves.
"""
struct PlaneWaveExpansion{P <: PropagationType, S <: SphereSamplingStrategy, C} <: AntennaFieldRepresentation{P, C}
   samplingstrategy :: S
   EθEϕ :: Matrix{C}
   wavenumber :: Number
end

"""
    PlaneWaveExpansion(Type{P}, samplingstrategy::S, Eθ::Matrix{C}, Eϕ::Matrix{C}, wavenumber::Number)

Returns a collection of electromagnetic plane waves propagating into various directions.

# Arguments:
- `Type{P}`: PropagationType
- `samplingstrategy::SphereSamplingStrategy`: sampling strategy which defines the propagation directions of the plane wave samples
- `Eθ::Matrix{C}`: θ-component amplitudes of the plane waves
- `Eϕ::Matrix{C}`: ϕ-component amplitudes of the plane waves
- `wavenumber`: wavenumber
"""
function PlaneWaveExpansion(P::PropagationType, samplingstrategy::S, Eθ::Matrix{C}, Eϕ::Matrix{C}, wavenumber::Number) where{S <: SphereSamplingStrategy, C}
   return PlaneWaveExpansion{P, S, C}(samplingstrategy, [Eθ Eϕ], wavenumber)
end
function _eθ(p::PlaneWaveExpansion)
   _ , s2 = _count_samples(p.samplingstrategy)
   return view(p.EθEϕ, :, 1:s2)
end
function _eϕ(p::PlaneWaveExpansion)
   _ , s2 = _count_samples(p.samplingstrategy)
   return view(p.EθEϕ, :, (s2 + 1) : 2*s2)
end
function asvector(p::PlaneWaveExpansion)
   return vec(p.EθEϕ)
end
Base.size(p::PlaneWaveExpansion) = size(asvector(p))
Base.getindex(p::PlaneWaveExpansion, i) = Base.getindex(asvector(p), i)
Base.setindex!(p::PlaneWaveExpansion, i, v) = Base.setindex(asvector(p), i, v)
function Base.similar(p::PlaneWaveExpansion{P, S, C}) where{P,S,C}
    return PlaneWaveExpansion{P, S, C}(p.samplingstrategy, similar(p.EθEϕ), E.wavenumber)
end
function setwavenumber!(p::PlaneWaveExpansion{P, S, C}, val) where{P,C,S}
    swe = PlaneWaveExpansion{P,C,S}(p.samplingstrategy, p.EθEϕ, val)
    return swe
end

function efield!(
   storage,
   pwe::PlaneWaveExpansion{Incident, S, C},
   R;
   reset=true
) where {S <: SphereSamplingStrategy, C<:Complex}
θweights, ϕweights, θs, ϕs= weightsandsamples(pwe.samplingstrategy)

T=real(C)
Rvec= SVector{3}(R)
reset && fill!(storage, zero(C))
   sintcost = sincos.(θs)
   for (kk, ϕ) in enumerate(ϕs)
       sinp, cosp = sincos(ϕ)
       eϕ = SVector{3}(-sinp, cosp, zero(T))
       for (k, θ) in enumerate(θs)
           sint, cost = sintcost[k]

           eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
           eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

           Epol = _eθ(pwe)[k, kk] * eθ + _eϕ(pwe)[k, kk] * eϕ
           ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
           storage += Epol * ejkr
       end
   end
   return storage
end

function hfield!(
   storage,
   pwe::PlaneWaveExpansion{Incident, S, C},
   R;
   reset=true
) where {S <: SphereSamplingStrategy, C<:Complex}
θweights, ϕweights, θs, ϕs= weightsandsamples(pwe.samplingstrategy)

T=real(C)
Rvec= SVector{3}(R)
reset && fill!(storage, zero(C))
   sintcost = sincos.(θs)
   for (kk, ϕ) in enumerate(ϕs)
       sinp, cosp = sincos(ϕ)
       eϕ = SVector{3}(-sinp, cosp, zero(T))
       for (k, θ) in enumerate(θs)
           sint, cost = sintcost[k]

           eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
           eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

           Hpol = _eθ(pwe)[k, kk] * eϕ - _eϕ(pwe)[k, kk] * eθ
           ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
           storage += Hpol * ejkr / Z₀
       end
   end
   return storage
end

function ehfield!(
   estorage,
   hstorage,
   pwe::PlaneWaveExpansion{Incident, S, C},
   R;
   reset=true
) where {S <: SphereSamplingStrategy, C<:Complex}
θweights, ϕweights, θs, ϕs= weightsandsamples(pwe.samplingstrategy)

T=real(C)
Rvec= SVector{3}(R)
reset && fill!(estorage, zero(C))
reset && fill!(hstorage, zero(C))
   sintcost = sincos.(θs)
   for (kk, ϕ) in enumerate(ϕs)
       sinp, cosp = sincos(ϕ)
       eϕ = SVector{3}(-sinp, cosp, zero(T))
       for (k, θ) in enumerate(θs)
           sint, cost = sintcost[k]

           eθ = SVector{3}(cosp * cost, sinp * cost, -sint)
           eᵣ = SVector{3}(cosp * sint, sinp * sint, cost)

           Hpol = _eθ(pwe)[k, kk] * eϕ - _eϕ(pwe)[k, kk] * eθ
           Epol = _eθ(pwe)[k, kk] * eθ + _eϕ(pwe)[k, kk] * eϕ
           ejkr = cis(-k₀ * udot(eᵣ, Rvec)) * θweights[k] * ϕweights[kk]
           storage .+= Epol * ejkr
           hstorage .+= Hpol * ejkr / Z₀
       end
   end
   return estorage, hstorage
end

function _add_or_reset!(
   storage::PlaneWaveExpansion,
   summand::PlaneWaveExpansion;
   reset::Bool = false,
)
   if reset
       storage.EθEϕ .= summand.EθEϕ
   else
      storage.EθEϕ .+= summand.EθEϕ
   end
   return storage
end
function _add_or_reset!(storage::AbstractMatrix, summand::AbstractMatrix; reset::Bool = false)
   if reset
       storage .= summand
   else
       storage .+= summand
   end
   return storage
end


function _muladd_or_mulreset!(
   storage::AbstractMatrix,
   summand::AbstractMatrix,
   factor::AbstractMatrix;
   reset::Bool = false,
)
   if reset
       storage .= factor .* summand
   else
       storage .= storage .+ (factor .* summand)
   end
   return storage
end
function _muladd_or_mulreset!(
   storage::PlaneWaveExpansion,
   summand::PlaneWaveExpansion,
   factor::Number;
   reset::Bool = false,
)
   if reset
       storage.EθEϕ .= factor .* summand.EθEϕ
   else
       storage.EθEϕ .+= factor .* summand.EθEϕ
   end
   return storage
end

function _mul_or_reset!(
   storage::PlaneWaveExpansion,
   factor::PlaneWaveExpansion;
   reset::Bool = false,
)
   if reset
       storage.EθEϕ .= factor.EθEϕ
   else
       storage.EθEϕ .*= factor.EθEϕ
   end
   return storage
end

function _mul_or_reset!(
   storage::PlaneWaveExpansion, factor::Number; reset::Bool = false)
   if reset
       fill!(storage.EθEϕ, factor)
   else
       storage.EθEϕ .*= factor
   end
   return storage
end
function _mul_or_reset!(factor::Number, storage::PlaneWaveExpansion; reset::Bool = false)
   _mul_or_reset!(storage, factor, reset = reset)
   return storage
end