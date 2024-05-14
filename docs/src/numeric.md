# Numeric types

This package exports a variety of types that represents elements of Clifford algebras.

## `AbstractCliffordNumber{Q,T}` and subtypes

The `AbstractCliffordNumber{Q,T}` type is the supertype for all implmentations of Clifford numbers.
`Q` is a `QuadraticForm`, which describes the number of dimensions with positive, negative, and zero
square, and `T` is a `Union{Real,Complex}` type of the coefficients.

!!! note "Future `StaticCliffordNumber{Q,T,L}` type"
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

These types represent Clifford numbers of exclusively even or odd grade, respectively. These are the
workhorses of geometric algebra, as they are produced through products of even or odd numbers of
1-blades. In the majority of cases, you can rely entirely on these types.

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
optimized so that operations with it are simply converted to scalar operations. Many operations on
Clifford numbers that return scalars will return a `KVector{0}` to preserve the metric signature and
other semantics associated with `AbstractCliffordNumber`.

!!! warning
    It is important to note that k-vectors are not *k-blades* (the wedge product of k 1-vectors) or
    *k-versors* (the geometric product of k 1-vectors). In dimensions up to 3, all k-vectors are
    also k-blades, but this is not generally true: as a counterexample, $e_1 e_2 + e_3 e_4$ is not
    representable as a k-blade. However, all k-blades are k-vectors. k-versors usually consist of
    more than one grade.

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
converting those to `CliffordNumber`. New `AbstractCliffordNumber` subtypes will widen directly to
`CliffordNumber` by default, and this should be overridden for new types so that it widens to the
smallest wider type.

## Construction and conversion

### From the constructors

The constructors of all `AbstractCliffordNumber` subtypes accept `Tuple` or `Vararg` arguments. In
interactive use, you will probably use the latter. When defining a new type, you only need to define
the `(::Type{T})(::Tuple)` constructors, as the `Vararg` constructors are automatically provided.

Some type parameters may be omitted in constructors, and the differences in behavior between these
constructors is given below, using `CliffordNumber` as an example:
- `CliffordNumber{Q,T}(x...)` converts all arguments `x` to type `T`.
- `CliffordNumber{Q}(x...)` promotes all arguments `x` to a common type `T`, so it is equivalent to `CliffordNumber{Q,promote_type(typeof.(x)...)}(x)`.
- `CliffordNumber(x...)` is *not* a valid constructor, as an algebra must be specified.

For types that include grade information, such as `KVector`, this information must be included to
produce a valid constructor.

#### Indices

In most literature, the components of multivectors are listed in grade order. However, this ordering
is not used here: instead, the natural binary ordering of blades is used.For a computer, each basis
blade of an ``n``-dimensional algebra can be represented with an ``n``-bit integer: each of its
binary digits correspond to the presence or absence of a vector composing the blade.

For a concrete example, the coefficients of a `CliffordNumber{VGA(3)}` are ordered like so:
```math
\left(1, e_1, e_2, e_1 e_2, e_3, e_1 e_3, e_2 e_3, e_1 e_2 e_3\right)
```
Note how ``e_1 e_2`` precedes ``e_3`` here, but also note how the first four elements and the second
four elements only differ by the presence of an `e_3` factor.

`CliffordNumber` and its backing `Tuple` can be indexed  straightforwardly with this relationship.
The basis blade order of all `AbstractCliffordNumber` instances are identical, with smaller types
like `KVector{2,VGA(3)}` skipping over all basis blades not of grade 2 in the list above.

!!! warning
    Many resources do not use a lexicographic order for the bivectors of the algebra of physical
    space or the spacetime algebra, opting for cyclic permutations so that ``e_3 e_1`` is preferred
    over ``e_1 e_3``. It's a good idea to check the convention before construction so you can
    include any necessary negative signs.

As a workaround for the possibly unintuitive ordering of coefficients, you can also use sums of
`KVector` instances: the sum will automatically promote to `EvenCliffordNumber`,
`OddCliffordNumber`, or `CliffordNumber` as needed to represent all grades.

#### Scalars

`CliffordNumber{Q}` and `EvenCliffordNumber{Q}` also accept a single scalar argument, and this
constructs an object with all non-scalar blade coefficients being zero. By definition, `KVector{0}`
does the same, and is the most efficient representation of a scalar associated with an algebra.

### From other Clifford numbers

Clifford numbers can be constructed from other CliffordNumbers. This implicitly performs a grade
projection operation, so this construction will always succeed, even if some of the basis blades of
the input are lost. By contrast, conversion will throw an `InexactError` if the result does not
contain all of the basis blades of the result.

```
julia> test = CliffordNumber{VGA(3)}(1, 2, 3, 4, 5, 6, 7, 8)
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
    than the input performs a grade projection operation without throwing an error.** However,
    *conversion* will throw an error if the grades of the input value are not present in the input
    type.
    
    **This is not how other subtypes of `Number` defined by Julia Base behave**, as their conversion
    operations are generally defined to be identical to the constructor, and always throw the same
    error for a given pair of type and value.

    If converting an `AbstractCliffordNumber` to any other numeric type, construction and conversion
    behave identically, as expected.

Construction and conversion of Clifford numbers from other Clifford numbers the only time that the
quadratic form type parameter can be omitted, as it can be inferred directly from the input. In the
case of `CliffordNumbers.Z2CliffordNumber`, the parity type parameter can also be inferred from a
`KVector` input.

### Scalar conversion

It may be desirable to convert the scalar type of a Clifford number without having to specify the
full typename of the desired output type. The function `scalar_convert(T, x)` takes a type 
`T<:Union{Real,Complex}` and any Clifford number `x` and converts its scalar type to
`T`. If `x` is a `Real` or `Complex`, it just converts `x` to an instance of `T`.
