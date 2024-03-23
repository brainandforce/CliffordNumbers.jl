# Numeric types

This package exports a variety of types that represents elements of Clifford algebras.

## `AbstractCliffordNumber{Q,T}` and subtypes

The `AbstractCliffordNumber{Q,T}` type is the supertype for all implmentations of Clifford numbers.
`Q` is a `QuadraticForm`, which describes the number of dimensions with positive, negative, and zero
square, and `T` is a `Union{Real,Complex}` type of the coefficients.

!!! note "Future `StaticCliffordNumber{Q,T,L}` type
    We may introduce a new abstract type, `StaticCliffordNumber{Q,T,L}`, for static implementations,
    like all of the ones provided by this package. These should be implemented as fixed length data
    structures (ideally an `NTuple{L,T}`).

### `CliffordNumber{Q,T,L}`: full grade Clifford numbers

`CliffordNumber{Q,T,L}` is the largest possible representation of a Clifford number, and it 
explicitly includes the coefficients for all `2^dimension(Q)` basis blades.

While this type is useful if working with objects that mix even and odd grades (for instance,
projectors or left minimal ideals), it is often more efficient to work with a smaller type, like the
ones described below.

### `EvenCliffordNumber{Q,T,L}` and `OddCliffordNumber{Q,T,L}`: even and odd graded elements

These types represent Clifford numbers of exclusively even or odd grade, respectively.

Internally, these are the same type: they alias `CliffordNumbers.Z2CliffordNumber{P,Q,T,L}`, where
`P` is a Boolean parameter which is `false` for `EvenCliffordNumber` and `true` for
`OddCliffordNumber`.

### `KVector{K,Q,T,L}`: elements of homogeneous grade

This type represents a k-vector, or a Clifford number of homogeneous grade, with the parameter `K`
indicating the grade.

It should be noted that in general, this type is not as efficient as `EvenCliffordNumber,`
`OddCliffordNumber`, or `CliffordNumber` when calculating products (though this may change in the 
future). Use this type if your primary operations are addition, or if you need compact storage.

However, there is one exception to this: `KVector{0}`, which represents a scalar. This type has been
optimized so that operations with it are simply converted to scalar operations.

!!! warning
    It is important to note that k-vectors are not *k-blades* (the wedge product of k 1-vectors) or
    *k-versors* (the geometric product of k 1-vectors). In dimensions up to 3, all k-vectors are
    also k-blades, but this is not generally true: as a counterexample, $e_1 e_2 + e_3 e_4$ is not
    representable as a k-blade. However, all k-blades are k-vectors.

!!! note
    In the future, we may consider adding a `DualKVector` or `PseudoKVector` type to more easily
    represent pseudoscalars, pseudovectors, and related objects.

## Promotion and widening

This package provides a robust promotion system for converting Clifford numbers and scalars to
common types before performing common numeric operations.

When promoting the types of Clifford numbers, there are two different types of promotions that can
occur: *scalar promotions*, which promote all the scalar types of the arguments to a common scalar
type, and *grade promotions*, which promote all the arguments to types which have a common set of
grades. This package provides the `scalar_promote` function that allows for the scalar types of each
argument to be promoted to a common type. Promote rules have been defined so that `promote` performs
a scalar promotion and a grade promotion. No function currently promotes only the grades of the
inputs.

The `widen` function in Julia Base widens an `Number` type to a type that can represent the result
of addition or subtraction with the the input type without overflowing or losing precision. This
functionality is passed through to Clifford numbers, but it only affects the scalar type, not the
grades.

The `widen_grade` function performs an equivalent operation with the grades of a Clifford number,
converting `KVector` to `EvenCliffordNumber` or `OddCliffordNumber` depending on the grade, and
converting those to `CliffordNumber`. New `AbstractCliffordNumber` subtypes should define this if
they intend to promote to types other than `CliffordNumber`.


## Construction and conversion

Clifford numbers can be constructed from other CliffordNumbers. This implicitly performs a grade
projection operation, so this construction will always succeed, even if some of the basis blades of
the input are lost. By contrast, conversion will throw an `InexactError` if the result does not
contain all of the basis blades of the result.

```
julia> test = CliffordNumber{APS}(1, 2, 3, 4, 5, 6, 7, 8)
8-element CliffordNumber{VGA(3), Int64}:
1 + 2σ₁ + 3σ₂ + 5σ₃ + 4σ₁σ₂ + 6σ₁σ₃ + 7σ₂σ₃ + 8σ₁σ₂σ₃

julia> EvenCliffordNumber(test)
4-element EvenCliffordNumber{VGA(3), Int64}:
1 + 4σ₁σ₂ + 6σ₁σ₃ + 7σ₂σ₃

julia> convert(EvenCliffordNumber, test)
ERROR: InexactError: ...
```

!!! danger
    This is an extremely important point: **construction of a Clifford number type with fewer grades
    than the input performs a grade projection operation.** Conversion *will* throw an error if the
    result is not exactly representable. 
    
    **This is not how other subtypes of `Number` defined by Julia Base behave**, as their conversion
    operations are generally defined to be identical to the constructor.

It may be desirable to convert the scalar type of an argument. The function `scalar_convert(T, x)`
takes a type `T<:Union{Real,Complex}` and any Clifford number `x` and converts its scalar type to
`T`. If `x` is a `Real` or `Complex`, it just converts `x` to an instance of `T`.
