# Indexing

Indexing is a critical operation in linear algebra: it would be difficult to imagine defining
operations on vectors, matrices, and arrays without some way of determining the coefficient at a
particular location in the object. Similarly, we often want to extract coefficients from
multivectors to perform operations.

However, the philosophy of this package -- treating multivectors a number system on the same footing
as that of complex numbers or quaternions -- means that the `AbstractCliffordNumber` type foregoes
array semantics. All `AbstractCliffordNumber` instances broadcast as scalars: for two instances `x`
and `y`, `x * y` is identical to `x .* y`, both of which calculate the geometric product.

On top of this, the grade representation of multivectors makes it difficult to relate the indices of
the backing `Tuple` for each type to the blades represented by the multivector for anything other
than the dense `CliffordNumber`. To solve this issue, this package provides types specifically
intended for indexing `AbstractCliffordNumber` subtypes: `BitIndex{Q}` and the subtypes of
`AbstractBitIndices{Q,C}`, `BitIndices{Q,C}` and `TransformedBitIndices{Q,C}`.

## `BitIndex{Q}`

The `BitIndex{Q}` type is the index type for `AbstractCliffordNumber{Q}`. This types wraps a `UInt`:
the first `dimension(Q)` bits represent the presence or absence of the basis vectors associated with
each dimension that are used to construct the indexed blade. The most significant bit is a sign bit:
this is needed to represent the parity associated with the order of wedge products used to construct
the blade.

Regardless of what position in the underlying `Tuple` a coefficient may be, if `x == y`, then the
same `BitIndex` `b` will index `x` and `y` identically:
```
julia> k = KVector{1,VGA(3)}(4,2,0)
3-element KVector{1, VGA(3), Int64}:
4e₁ + 2e₂

julia> b = BitIndex(k, 1)
BitIndex(Val(VGA(3)), 1)

julia> k == OddCliffordNumber(k)
true

julia> k[b] == OddCliffordNumber(k)[b]
true
```
It should also be noted that indexing an `AbstractCliffordNumber` at an index which is not
explicitly represented by the type returns zero:
```
julia> k[BitIndex(k, 1, 2)]
0

```

### Construction

The internal constructor `BitIndex{Q}(signbit::Bool, blade::Unsigned)` converts `blade` to a `UInt`
and changes the most significant bit to match `signbit`. However, in many cases, constructing a
`BitIndex{Q}` directly from the sign bit and representation of the blade as an unsigned integer is
inconvenient.

For this reason, we define two constructors, the first being `BitIndex(::Val{Q}, i::Integer...)`,
which takes the algebra `Q` wrapped in `Val` (for reasons of type stability) and the integers `i`
corresponding to basis blade indices defined in `Q`. The second constructor, 
`BitIndex(x, i::Integer...)` calls `BitIndex(Val(signature(x))), i...)`, automatically determining
`Q`.

In both cases, the parity of the permutation of indices is determined automatically from the integer
arguments, so this bit is automatically assigned. The presence of a parity bit allows for correct 
indexing with non-lexicographic conventions: in some literature, ``e_3 e_1`` is used instead of
``e_1 e_3`` for one of the bivector components of APS and STA.
```
julia> l = KVector{2,VGA(3)}(0,6,9)
3-element KVector{2, VGA(3), Int64}:
6e₁e₃ + 9e₂e₃

julia> l[BitIndex(l, 1, 3)]
6

julia> l[BitIndex(l, 3, 1)]
-6
```

### Tuples of `BitIndex{Q}`

Julia `AbstractArray` instances can be indexed not just with integers, but with arrays of integers
or special objects representing iterable ranges, such as `:`. Perhaps surprisingly, tuples 
containing integers (or any other valid index object) are not valid indices of `AbstractArray`:
```
julia> (1:10)[(2, 3, 4)]
ERROR: ArgumentError: invalid index: (2, 3, 4) of type Tuple{Int64, Int64, Int64}

```

In the case of `AbstractCliffordNumber`, we have a compelling reason to use tuples of `BitIndex{Q}`
objects for indexing: since the length of a `Tuple` is known statically, we can use that information
to construct a new `Tuple` of coefficients with statically known length, which may be useful if we
want to leverage indexing to convert types. Therefore, indexing an `AbstractCliffordNumber{Q,T}` 
with an `NTuple{N,BitIndex{Q}}` returns an `NTuple{N,T}`, which can be fed into the constructor for 
a different type, and this is what the package uses internally to perform conversion.

### Operations

`BitIndex{Q}` supports a variety of unary and binary operations, many of which are used internally
for tasks like calculating geometric products. Many of these operations are also supported for
`AbstractCliffordNumber{Q}` instances, such as negation (`-`) and the geometric product (`*`).

## `AbstractBitIndices{Q,C}`

Considering that the indices of an `AbstractCliffordNumber` provided by this package are known from
the type, it makes sense to define a type which represents all indices of a subtype or instance of
`AbstractCliffordNumber`. This package defines the supertype `AbstractBitIndices{Q,C}` and concrete
types `BitIndices{Q,C}` and `TransformedBitIndices{Q,C,F}` that allow for easy operation across all
valid indices of a type `C`.

### `BitIndices{Q,C}`

The `BitIndices` object for an `AbstractCliffordNumber` subtype `C` or instance `x` can be 
constructed with `BitIndices(C)` or `BitIndices(x)`. This is a singleton type which is also returned
by `eachindex(x)`.

The type information of the `BitIndices` object is used to determine the type of the result of 
indexing. Conveniently, you can index one Clifford number `x` with the indices of another Clifford 
number `y`, and this converts `x` to the a type similar to that of `y`, but retaining the element
type of `x`:
```
julia> k = KVector{1,VGA(3)}(4,2,0)             # eltype Int64
3-element KVector{1, VGA(3), Int64}:
4e₁ + 2e₂

julia> z = zero(CliffordNumber{VGA(3),Float64}) # eltype Float64
8-element CliffordNumber{VGA(3), Float64}:
0.0

julia> k[BitIndices(ans)]
8-element CliffordNumber{VGA(3), Int64}:
4e₁ + 2e₂

```

`BitIndices{Q,C}` subtypes `AbstractBitIndices{Q,C}`, which subtypes `AbstractVector{BitIndex{Q}}`,
meaning that the objects can be indexed with integers. This indexing is one-based to match the
indices of the `Tuple` backing a `CliffordNumber`, but this is subject to change in the future.

!!! warning "Avoiding type proliferation with `BitIndices`"
    One fundamental issue with `BitIndices{Q,C}` is that `C` is only constrained to be a subtype of
    `AbstractCliffordNumber{Q}`. However, `C` may have other type parameters that vary, but produce
    objects that index equivalently and are identical elementwise:
    ```
    julia> BitIndices{STA,CliffordNumber{STA}}() == BitIndices{STA,CliffordNumber{STA,Float32}}()
    true

    julia> BitIndices{STA,CliffordNumber{STA}}() == BitIndices{STA,CliffordNumber{STA,Float32,8}}()
    true 
    ```
    For this reason, you should call `BitIndices(x)` or `BitIndices(C)` instead of directly
    constructing `BitIndices{Q,C}()`.

### `TransformedBitIndices{Q,C,F}`

If we want to perform a transformation on all elements of a `BitIndices` instance, we may use dot 
syntax and broadcast an operation: for instance, `reverse.(BitIndices(x))` to obtain the reverse of
all indices. Without any explicit specification of behavior, this will return a `Vector`, and for 
the sake of performance we'd like to avoid returning types without statically known lengths.

We could solve this by converting `BitIndices` to a `Tuple` internally as part of the broadcasting
implementation, but doing this will cause its own set of problems: indexing an
`AbstractCliffordNumber` with a `Tuple` returns a `Tuple`. This may be fine internally, if we know
to call the constructor, but this is potentially confusing for new users. For this reason, we 
provide `TransformedBitIndices{Q,C,F}`, which is a lazy representation of a function applied to
every element of `BitIndices(C)`.

For convenience, we provide a few aliases for operations which are commonly used with `BitIndices`:
- `ReversedBitIndices{Q,C}` is the type of `reverse.(::BitIndices{Q,C})`.
- `GradeInvolutedBitIndices{Q,C}` is the type of `grade_involution.(::BitIndices{Q,C})`.
- `ConjugatedBitIndices{Q,C}` is the type of `conj.(::BitIndices{Q,C})`.

In the future, we may override some `broadcast` implementations to ensure that all of these types
are interconvertible with each other and with `BitIndices{Q,C}`.
