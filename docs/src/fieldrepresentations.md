# [Representation of Antenna Fields](@id fieldrepresentation)

The electromagnetic fields of an antenna can be represented in various ways. Some of the most important field representations known from literature are implemented in `AntennaFieldRepresentations.jl`. 

The core of `AntennaFieldRepresentations.jl` revolves around implemantations of the abstract type `AntennaFieldRepresentation`. Implementations of the type [`AntennaFieldRepresentation`](@ref) provide equivalent representations of the electromagnetic fields of an antenna (in a certain region of space where the representation converges). Each instance of `AntennaFieldRepresentation` is a discretized version of the antenna fields, i.e., a collection of coefficients from which the antenna fields can be (approximately) calculated. One of the imperatives behind `AntennaFieldRepresentations.jl` is that an `AntennaFieldRepresentation` behaves like an `AbstractVector{C}` with extra context. In particular, the wavenumber of the radiation frequency of the antenna field is stored in the `AntennaFieldRepresentation`.

!!! tip
    Implementations of the type `AntennaFieldrepresentation` support the interface of `Base.AbstractVector`. In particular the methods `Base.getindex` and `Base.setindex!` are supported which allow to access the coefficient vector by indexing the struct representing the fields directly by square brackets `[]`.
---

## Three Different Propagation Types of `AntennaFieldRepresentation`

`AntennaFieldRepresentation`s come in up to three different propagation types.
For indicating the propagation type of a concrete `AntennaFieldRepresentation`-type, we use the singleton-types 
- `Radiated`,
- `Absorbed`,
- or `Incident`,
although not every implementation of `AntennaFieldRepresentation` will support all three variants. 

Internally, `Radiated`, `Absorbed`, or `Incident` are represented by singleton types, which all subtype the abstract type [`PropagationType`](@ref). The `PropagationType` is stored as the parameter `P` in the type definition of `AntennaFieldRepresentation{P, C} ` to allow for efficient dispatch of methods implementing the `AntennaFieldRepresentation` interface (the other parameter `C` denotes the element type of the coefficient vector).

## `AntennaFieldRepresentation` Interface

Any subtype of `AntennaFieldrepresentation` implements the following methods[^1]:

| Method name               | Optional | Fallback method      | Short Description                                                     |
| :------------------------ |:-------- | :------------------- | :-------------------------------------------------------------------- |
| `Base.size`               | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.getindex`           | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.setindex!`          | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.similar`            | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| [`asvector`](@ref)        | No       |                      | Return the coefficients as a vector                                   |
|                           |          |                      |                                                                       |
| [`efield`](@ref)          | Yes      | `efield!`            | Return the E-field vector at given location                           |
| [`efield!`](@ref)         | No       |                      | In-place version of `efield`                                          |
| [`hfield`](@ref)          | Yes      | `hfield!`            | Return the H-field vector at given location                           |
| [`hfield!`](@ref)         | No       |                      | In-place version of `hfield`                                          |
| [`ehfield`](@ref)         | Yes      | `ehfield!`           | Simulataneous calculation of E-field and H-field                      |
| [`ehfield!`](@ref)        | Yes      | `efield!`, `hfield!` | In-place version of `ehfield`                                         |
| [`farfield`](@ref)        | No       |                      | Return the farfield vector at given direction                         |
|                           |          |                      |                                                                       |
| [`getwavenumber`](@ref)   | No       |                      | Return the wavenumber of the radiation frequency                      |
| [`setwavenumber!`](@ref)  | No       |                      | Set the wavenumber of the radiation frequency                         |
| [`rotate`](@ref)          | Yes      | `rotate!`            | Rotate the representation by the Euler angles `χ`, `θ`, `ϕ`           |
| [`rotate!`](@ref)         | No       |                      | In-place version of `rotate`                                          |
| [`spatialshift`](@ref)    | Yes      | `spatialshift!`      | Spatially shift representation to new location `R`                    |
| [`spatialshift!`](@ref)   | No       |                      | In-place version of `spatialshift`                                    |
|                           |          |                      |                                                                       |
| [`equivalentorder`](@ref) | No       |                      | Return the estimated spherical mode order L of the representation     |
|                           |          |                      |                                                                       |
| [`changerepresentation`](@ref)    | Yes      |              | Change the representation type. Not all resulting types are supported |
| `transmit`                | No       |                      | Calculate transmission between AUT and field sampling                 |

[^1]: Note for developers: If you want to implement your own subtype of `AntennaFieldrepresentation`, make sure to support this interface to adhere to the general functionality of this package.
