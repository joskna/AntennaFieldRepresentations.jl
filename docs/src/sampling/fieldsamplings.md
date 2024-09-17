# [Sampling Setups for Antenna Fields](@id fieldsampling)

In electromagnetic field solvers a concept called *field monitor* or *field request* (or similiar) is very common: The desired electromagnetic fields are not calculated everywhere in space but only on a limited set of points defined in the *field request*. This saves resources and makes the computation of the electromagnetic fields feasible in the first place.

In `AntennaFieldRepresentations.jl`, a similar concept exists. A `FieldSampling` is a struct which defines the positions in which the electromagnetic fields of an [`AntennaFieldrepresentation`](@ref fieldrepresentation) shall be sampled. Furthermore, many `Fieldsampling`s allow to take the receive pattern of a realistic probe antenna into account. This way, rather complex and realistic antenna measurement setups can be modelled.

A `FieldSampling` can be used to create a [`TransmissionMap`](@ref transmissionmap), a functional relationship between the coefficients of an [`AntennaFieldRepresentation`](@ref fieldrepresentation) and the output signals of the probe antennas defined by the `FieldSampling`[^1]. 
The [`TransmissionMap`](@ref transmissionmap) serves as a linear operator (a matrix) between the coefficient vector of the [`AntennaFieldRepresentation`](@ref fieldrepresentation) and the vector of output signals of the probe antennas.
Depending on the specific types of the [`AntennaFieldRepresentation`](@ref fieldrepresentation) and the `FieldSampling`, the algorithm implementing the matrix-vector product for the `TransmissionMap` will change. 
Particularly efficient algorithms arise for certain combinations of [`AntennaFieldRepresentation`](@ref fieldrepresentation)s and `FieldSampling`s. 

Every `FieldSampling` reserves a storage for the measured values of the ``S_{21}``-parameter between the probe antennas and an [`AntennaFieldRepresentation`](@ref fieldrepresentation). Internally, the storage can have various forms to best match the situation described by the `FieldSampling` implementation. However, all `FieldSampling`s implement a common interface which let the user access the storage for the measured ``S_{21}``-parameters a [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).

!!! tip
    Implementations of the type `FieldSampling` support the interface of `Base.AbstractVector`. In particular the methods `Base.getindex` and `Base.setindex!` are supported which allow to access the measurement vector by indexing the struct representing the fields directly by square brackets `[]`.
---

## `FieldSampling` Interface

Any subtype of `FieldSampling` implements the following methods[^2]:

| Method name               | Optional | Fallback method      | Short Description                                                     |
| :------------------------ |:-------- | :------------------- | :-------------------------------------------------------------------- |
| `Base.size`               | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.getindex!`          | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.setindex`           | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `Base.similar`            | No       |                      | See interface of [`Base.AbstractVector`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)|
| `asvector`                | No       |                      | Return the measurement samples as a vector                            |


[^1]: Hertzian dipoles (short electric dipoles) and Fitzgerald dipoles (short magnetic dipoles) can formally be used as probe antennas to directly sample electric and magnetic fields of an `AntennaFieldRepresentation`: Their output signals are directly proportional to the electric and magnetic field component which is parallel to their dipole orientation, respectively.

[^2]: Note for developers: If you want to implement your own subtype of `FieldSampling`, make sure to support this interface to adhere to the general functionality of this package. 