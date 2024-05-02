# Metric signatures

Clifford algebras are characterized by the metric signatures of the spaces they represent. Some are
very commonly used, such as the algebra of physical space (APS), or are generated as part of a
family, such as the projective geometric algebras (PGAs), but in other cases you may need the 
flexibility to work with custom metric signatures.

The `CliffordNumbers.Metrics` submodule provides tools for working with metric signatures.

## Interface

The type parameter `Q` of `AbstractCliffordNumber{Q,T}` is not constrained in any way, which means
that any type or data consisting of pure bits may reside there. However, for the sake of
correctness and fully defined behavior, `Q` must satisfy an informal interface.

Metric signature objects are treated like `AbstractVector{Int8}` instances, but with the elements
constrained to be equal to `+1`, `0`, or `-1`, corresponding to basis 1-blades squaring to positive
values, negative values, or zero. In the future, we may support arbitrary values for this type.

This array is not constrained to be a 1-based array, and the values of `eachindex` for the array
correspond to the indices of the basis 1-blades of the algebra.

## The `Metrics.AbstractSignature` type

We define a type, `Metrics.AbstractSignature <: AbstractVector{Int8}`, for which this interface is
already partially implemented.

## Pre-defined signatures

There are many commonly used families of algebras, and for the sake of convenience, we provide four
subtypes of `Metrics.AbstractSignature` to handles these cases:

  * `Metrics.VGA` represents vanilla geometric algebras.
  * `Metrics.PGA` represents projective geometric algebras.
  * `Metrics.CGA` represents conformal geometric algebras.
  * `Metrics.LGA{C}` represents Lorentzian geometric algebras:
      * `Metrics.LGAEast` uses the East Coast convention (timelike dimensions square to -1).
      * `Metrics.LGAWest` uses the West Coast convention (timelike dimensions square to +1).

To construct an instance of one of these types, call it with the number of modeled spatial
dimensions:
  * `Metrics.VGA(3)` models 3 spatial dimensions with no extra dimensions.
  * `Metrics.PGA(3)` models 3 spatial dimensions with 1 degenerate (zero-squaring) dimension.
  * `Metrics.CGA(3)` models 3 spatial dimensions with 2 extra dimensions.
  * `Metrics.LGAEast(3)` models 3 spatial dimensions with an extra negative-squaring time dimension.
