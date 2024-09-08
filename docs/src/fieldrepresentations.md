# Representation of Antenna Fields

The electromagnetic fields of an antenna can be represented in various ways. Some of the most important field representations known from literature are implemented in `AntennaFieldRepresentations.jl`. 

The core of `AntennaFieldRepresentations.jl` revolves around implemantations of the abstract type `AntennaFieldRepresentation`. Implementations of the type `AntennaFieldRepresentation` provide equivalent representations of the electromagnetic fields of an antenna (in a certain region of space where the representation converges). Each instance of `AntennaFieldRepresentation` is a discretized version of the antenna fields, i.e., a collection of coefficients from which the antenna fields can be (approximately) calculated. One of the imperatives behind `AntennaFieldRepresentations.jl` is that an `AntennaFieldRepresentation` behaves like an `AbstractVector{C}` with extra context.

!!! tip
    Implementations of the type `AntennaFieldrepresentation` support the interface of `Base.AbstractVector`. In particular the methods `Base.getindex` and `Base.setindex!` are supported which allow to access the coefficient vector by indexing the struct representing the fields directly by square brackets `[]`.
---

Any subtype of `AntennaFieldrepresentation` implements the following methods[^1]:

| Method name               | Optional | Fallback method      | Short Description                                                     |
| :------------------------ |:-------: | :------------------- | :-------------------------------------------------------------------- |
| `Base.size`               | No       |                      | See interface of `Base.AbstractVector`                                |
| `Base.getindex!`          | No       |                      | See interface of `Base.AbstractVector`                                |
| `Base.setindex`           | No       |                      | See interface of `Base.AbstractVector`                                |
| `Base.similar`            | No       |                      | See docimentations of `Base.similar`                                  |
| `as_vector`               | No       |                      | Return the coefficient vector                                         |
|                           |          |                      |                                                                       |
| `efield`                  | Yes      | `efield!`            | Return the E-field vector at given location                           |
| `efield!`                 | No       |                      | In-place version of `efield`                                          |
| `hfield`                  | Yes      | `hfield!`            | Return the H-field vector at given location                           |
| `hfield!`                 | No       |                      | In-place version of `hfield`                                          |
| `ehfield`                 | Yes      | `ehfield!`           | Simulataneous calculation of E-field and H-field                      |
| `ehfield!`                | Yes      | `efield!`, `hfield!` | In-place version of `ehfield`                                         |
| `farfield`                | No       |                      | Return the farfield vector at given direction                         |
|                           |          |                      |                                                                       |
| `rotate`                  | Yes      | `rotate!`            | Rotate the representation by the Euler angles `χ`, `θ`, `ϕ`           |
| `rotate!`                 | No       |                      | In-place version of `rotate`                                          |
| `spatialshift`            | Yes      | `spatialshift!`      | Spatially shift representation to new location `R`                    |
| `spatialshift!`           | No       |                      | In-place version of `spatialshift`                                    |
|                           |          |                      |                                                                       |
| `changerepresentation`    | Yes      |                      | Change the representation type. Not all resulting types are supported |
| `transmit`                | No       |                      | Calculate transmission between AUT and field sampling                 |

[^1]: Note for developers: If you want to implement your own subtype of `AntennaFieldrepresentation`, make sure to support this interface to adhere to the general functionality of this package.