# Operation Maps for Operations on Antenna Field Representations
In `AntennaFieldRepresentations.jl`, an operation map is a [function-like object](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects-1) which corresponds to a certain method.
All operation maps  in `AntennaFieldRepresentations.jl` subtype the abstract type `OperationMap{C}`. 
The idea behind these `OperationMap`s is to reuse allocated memory and save redundant operations for repeated method calls.
The workflow involving `OperationMap`s starts with initializeing an instance of the `OperationMap`. 
The generated object then can be used like a function, where all overhead for the function call is moved to the initialization of the object. 

!!! todo
    Example coming soon ...
---

For many methods in `AntennaFieldRepresentations.jl`, a correspoding operation map is defined. 
A comprehensive list of all operation maps provided by `AntennaFieldRepresentations.jl` and their corresponding methods is given in the table below.

| Operation Map                 | Corresponding method      | Short Description                                         |
| :---------------------------- | :------------------------ | :-------------------------------------------------------- |
| `TransmitMap`                 | `transmit`                | Transmission operator between an `AntennaFieldRepresentation` and a `Fieldsampling` |
| `RotateMap`                   | `rotate`                  | Rotation operator for a certain `AntennaFieldRepresentation` around fixed Euler angles ``\vartheta``, ``\varphi``,  and ``\chi``|
| `SpatialShiftMap`             | `spatialshift`            | Shift operator for a certain `AntennaFieldRepresentation` which moves the `AntennaFieldRepresentation` along a vector `R` |
| `TransferMap`                 | `transfer`                | Transfer operator which transfers a certain `AntennaFieldRepresentation` of `Radiated` type into a corresponding `Incident` representation in a different coordinate system which is shifted by a vector `R` from the original coordinate system |
| `ChangeRepresentationMap`     | `changerepresentation`    | Conversion operator from one type of `AntennaFieldrepresentation` into another|
| `ResampleMap`              | `resample`             | Resample operator of a `PlaneWaveExpansion` into a `PlaneWaveExpansion` with a different `SphereSamplingStrategy`|

## [Operation Maps as Linear Operators](@id operationmaps_linmap)
The abstract type `OperationMap{C} <: LinearMaps.LinearMap{C}` is a subtype of [`LinearMaps.LinearMap`](https://julialinearalgebra.github.io/LinearMaps.jl/stable/generated/custom/) which means that the operations represented by an `OperationMap{C}` are strictly linear operations.
The usual way to evaluate an `OperationMap` is to use it in a matrix-vector-product with an `AntennaFieldRepresentation`. 

!!! warning
    Calling an `OperationMap` with a different type of `AntennaFieldRepresentation` could lead to unexpected behavior or errors.
---

!!! todo
    Example coming soon ...
---


Since an `OperationMap{C}` implements the interface of [`LinearMaps.jl`](https://julialinearalgebra.github.io/LinearMaps.jl/stable/), the following methods are inherently supported with any operation map `M`[^1]

| Method                    | Optional | Short Description                                                  |
| :------------------------ |:-------- | :----------------------------------------------------------------- |
| `Base.size(M)`            | No       | Dimensions of the operator                                         |
| `M*x`                     | Yes      | Out of place multiplication                                        |
| `mul!(y, M, x)`           | Yes      | in-place multiplication with vectors                               |
| `mul!(y, M, x, α, β)`     | Yes      | in-place multiply-and-add with vectors                             |
| `mul!(Y, M, X, α, β)`     | Yes      | in-place multiplication and multiply-and-add with matrices         |
| `Matrix(M)`               | Yes      | Conversion to a matrix                                             |
| `M[:, i]`                 | Yes      | Complete slicing of columns (and rows if the adjoint action is defined) |
 
Remember that any `AntennaFieldRepresentation` is interpreted as an `AbstractVector` with additional context in `AntennaFieldRepresentations.jl`. Thus, by using the interface of [`LinearMaps.jl`](https://julialinearalgebra.github.io/LinearMaps.jl/stable/), the user inherently implies that the vectors `x` which were used in the above matrix vector products corresponds to the coefficient vector belonging to the `AntennaFieldRepresentation` used for the definition of the `OperationMap{C}`.

## Adjoint, Transposed, and Inverse Operation Maps
In `AntennaFieldRepresentations.jl`, it is easy to define the adjoint, transposed, or the inverse of any operation map. Generically, taking the transpose of an operation map wraps the operation map by a `TransposeMap`, taking the adjoint wraps it by an `AdjointMap`, and taking the inverse wraps it by an `InverseMap` (or an `IterativeInverseMap`). 


In particular the inverse operations can be expensive to compute directly. For many operation maps, the function call `inverse(M)` will result in an `IterativeInverseMap`, where the matrix vector product with the inverse operator is approximated with the help of [`IterativeSolvers.jl`](https://iterativesolvers.julialinearalgebra.org/stable/).


[^1]: Note for developers: If you want to implement your own operation map `M`, it is sufficient to implement the methods  `Base.size(M)` and `LinearMaps._unsafe_mul!(y, A::MyFillMap, x::AbstractVector)` for supporting the basic features.

