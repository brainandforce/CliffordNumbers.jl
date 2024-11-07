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
intended for indexing `AbstractCliffordNumber` subtypes: `BladeIndex{Q}` and `BladeIndices{Q,C}`.

## `BladeIndex{Q}`

The `BladeIndex{Q}` type is the index type for `AbstractCliffordNumber{Q}`. This types wraps a 
`UInt`: the first `dimension(Q)` bits represent the presence or absence of the basis vectors
associated with each dimension that are used to construct the indexed blade. The most significant
bit is a sign bit: this is needed to represent the parity associated with the order of wedge
products used to construct the blade.

Regardless of what position in the underlying `Tuple` a coefficient may be, if `x == y`, then the
same `BladeIndex` `b` will index `x` and `y` identically:
```julia-repl
julia> k = KVector{1,VGA(3)}(4,2,0)
3-element KVector{1, VGA(3), Int64}:
4e₁ + 2e₂

julia> b = BladeIndex(k, 1)
BladeIndex(Val(VGA(3)), 1)

julia> k == OddCliffordNumber(k)
true

julia> k[b] == OddCliffordNumber(k)[b]
true
```
It should also be noted that indexing an `AbstractCliffordNumber` at an index which is not
explicitly represented by the type returns zero:
```julia-repl
julia> k[BladeIndex(k, 1, 2)]
0
```
The only way of throwing an error with a `BladeIndex` object is by having a mismatch of algebra type
parameters between the Clifford number and the index.

### Construction

The internal constructor `BladeIndex{Q}(signbit::Bool, blade::Unsigned)` converts `blade` to a
`UInt` and changes the most significant bit to match `signbit`. However, in many cases, constructing
a `BladeIndex{Q}` directly from the sign bit and representation of the blade as an unsigned integer 
is inconvenient.

For this reason, we define two constructors, the first being `BladeIndex(::Val{Q}, i::Integer...)`,
which takes the algebra `Q` wrapped in `Val` (for reasons of type stability) and the integers `i`
corresponding to basis blade indices defined in `Q`. The second constructor, 
`BladeIndex(x, i::Integer...)` calls `BladeIndex(Val(signature(x))), i...)`, automatically
determining `Q`.

In both cases, the parity of the permutation of indices is determined automatically from the integer
arguments, so this bit is automatically assigned. The presence of a parity bit allows for correct 
indexing with non-lexicographic conventions: in some literature, ``e_3 e_1`` is used instead of
``e_1 e_3`` for one of the bivector components of APS and STA.
```
julia> l = KVector{2,VGA(3)}(0,6,9)
3-element KVector{2, VGA(3), Int64}:
6e₁e₃ + 9e₂e₃

julia> l[BladeIndex(l, 1, 3)]
6

julia> l[BladeIndex(l, 3, 1)]
-6
```

### Tuples of `BladeIndex{Q}`

Julia `AbstractArray` instances can be indexed not just with integers, but with arrays of integers
or special objects representing iterable ranges, such as `:`. Perhaps surprisingly, tuples 
containing integers (or any other valid index object) are not valid indices of `AbstractArray`, even
though vectors of integers are:
```julia-repl
julia> (1:10)[[2,3,4]]
3-element Vector{Int64}:
 2
 3
 4

julia> (1:10)[(2, 3, 4)]
ERROR: ArgumentError: invalid index: (2, 3, 4) of type Tuple{Int64, Int64, Int64}
```
In the case of `AbstractCliffordNumber`, we have a compelling reason to use tuples of
`BladeIndex{Q}` objects for indexing: since the length of a `Tuple` is known statically, we can use
that information to construct a new `Tuple` of coefficients with statically known length, which may
be useful if we want to leverage indexing to convert types. Therefore, indexing an
`AbstractCliffordNumber{Q,T}` with an `NTuple{N,BladeIndex{Q}}` returns an `NTuple{N,T}`, which can
be fed into the constructor for a different type, and this is what the package uses internally to
perform conversion.

### Operations

`BladeIndex{Q}` supports a variety of unary and binary operations, many of which are used internally
for tasks like calculating geometric products. Many of these operations are also supported for
`AbstractCliffordNumber{Q}` instances, such as negation (`-`) and the geometric product (`*`).

## `BladeIndices{Q,C}` and related types

Considering that the indices of an `AbstractCliffordNumber` provided by this package are known from
the type, it makes sense to define a type which represents all indices of a subtype or instance of
`AbstractCliffordNumber`. This package defines the singleton type `BladeIndices{Q,C}`, a subtype of
`AbstractArray{BladeIndex{Q}}`, which compactly represents the unique, explicitly present indices of
a Clifford number of type `C`.

`BladeIndices` objects can be constructed using `BladeIndices`, `eachindex`, or `keys`:
```julia-repl
julia> k = KVector{1,VGA(3)}(4, 2, 0)
3-element KVector{1, VGA(3), Int64}:
4e₁ + 2e₂

julia> BladeIndices(k)
3-element BladeIndices{VGA(3), KVector{1, VGA(3)}}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)

julia> eachindex(k)
3-element BladeIndices{VGA(3), KVector{1, VGA(3)}}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)

julia> keys(k)
3-element BladeIndices{VGA(3), KVector{1, VGA(3)}}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)
```
Notice how the output type lacks the scalar type parameter or length parameter of the type of the
`AbstractCliffordNumber`. Stripping the scalar type prevents excessive proliferation of indexing
types (which is useful for generated functions), and ensures that the indexing types of two Clifford
numbers with different scalar types are equal in the `===` sense:
```julia-repl
julia> kk = float(k)
3-element KVector{1, VGA(3), Float64}:
4.0e₁ + 2.0e₂

julia> BladeIndices(kk) === BladeIndices(k)
true
```
For this reason (and others), we recommend using `BladeIndices(k)`, `eachindex(k)`, or `keys(k)` as
appropriate in generic code instead of `BladeIndices{signature(k),typeof(k)}()`.

### Conversion

One of the greatest uses for `BladeIndices` objects is that they can be used to perform tasks such
as conversion or grade projection:
``` julia-repl
julia> inds = BladeIndices(OddCliffordNumber{VGA(3)})
4-element BladeIndices{VGA(3), OddCliffordNumber{VGA(3)}}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)
 BladeIndex(Val(VGA(3)), 1, 2, 3)

julia> k[inds]
4-element OddCliffordNumber{VGA(3), Int64}:
4e₁ + 2e₂
```
This reindexing operation converted `k` from a `KVector{1,VGA(3)}` to an `OddCliffordNumber{VGA(3)}`
automatically. All construction and conversion operations between Clifford numbers of the same 
algebra have a simple implementation in terms of indexing.

### Applying involutions

`BladeIndices{Q,C}` is actually not an independent type: it is an alias of an internal type,
`CliffordNumbers.CGNBladeIndices{Q,C,N,I,R}`. The final three type parameters lazily represent three
common operations: negation (`N`), grade involution (`I`), and reversion (`R`). Their composition is
efficiently represented in these type parameters.

Rather than deal directly with the possible types and their instances, custom broadcasting
implementations are used (though the return type is visible):
```julia-repl
julia> grade_involution.(inds)
4-element CliffordNumbers.CGNBladeIndices{VGA(3), OddCliffordNumber{VGA(3)}, false, true, false}:
 -BladeIndex(Val(VGA(3)), 1)
 -BladeIndex(Val(VGA(3)), 2)
 -BladeIndex(Val(VGA(3)), 3)
 -BladeIndex(Val(VGA(3)), 1, 2, 3)

julia> reverse.(inds) # equivalent to adjoint.(inds)
4-element CliffordNumbers.CGNBladeIndices{VGA(3), OddCliffordNumber{VGA(3)}, false, false, true}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)
 -BladeIndex(Val(VGA(3)), 1, 2, 3)
```
The Clifford conjugate is equivalent to combining grade involution with reversion:
```
julia> conj.(inds)
4-element CliffordNumbers.CGNBladeIndices{VGA(3), OddCliffordNumber{VGA(3)}, false, true, true}:
 -BladeIndex(Val(VGA(3)), 1)
 -BladeIndex(Val(VGA(3)), 2)
 -BladeIndex(Val(VGA(3)), 3)
 BladeIndex(Val(VGA(3)), 1, 2, 3)
```
And as with all involutions, applying them twice undoes them, leaving behind an ordinary
`BladeIndices` object:
```
julia> conj.(conj.(inds))
4-element BladeIndices{VGA(3), OddCliffordNumber{VGA(3)}}:
 BladeIndex(Val(VGA(3)), 1)
 BladeIndex(Val(VGA(3)), 2)
 BladeIndex(Val(VGA(3)), 3)
 BladeIndex(Val(VGA(3)), 1, 2, 3)
```

!!! warn
    `reverse(::BladeIndices)` actually reverses the elements of a `BladeIndices` object, just as it
    would with any other `AbstractVector`! The dot syntax for broadcasting is necessary:
    ```julia-repl
    julia> reverse(inds)
    4-element Vector{BladeIndex{VGA(3)}}:
     BladeIndex(Val(VGA(3)), 1, 2, 3)
     BladeIndex(Val(VGA(3)), 3)
     BladeIndex(Val(VGA(3)), 2)
     BladeIndex(Val(VGA(3)), 1)
    ```
    `adjoint(::BladeIndices)` does broadcast the reverse operation (it is equivalent to `'`), but
    also changes the array type to a `LinearAlgebra.Adjoint` wrapper:
    ```julia-repl
    julia> adjoint(inds)
    1×4 adjoint(::BladeIndices{VGA(3), OddCliffordNumber{VGA(3)}}) with eltype BladeIndex{VGA(3)}:
     BladeIndex(Val(VGA(3)), 1)  BladeIndex(Val(VGA(3)), 2)  BladeIndex(Val(VGA(3)), 3)  -BladeIndex(Val(VGA(3)), 1, 2, 3)
    ```
    
As with conversion and construction, efficient implementations of the common involutions are all
done in terms of indexing with `CGNBladeIndex` objects.
