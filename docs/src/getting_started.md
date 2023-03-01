# Getting started

## The `CliffordNumber` data type

A `CliffordNumber{Q,T,L}` is a Clifford number associated with a `QuadraticForm` `Q`, a backing
`Real` or `Complex` type `T`, and a length `L`. The length parameter is redundant, and in many
cases, it may be omitted without consequence.

!!! warning Although the length may be omitted in many cases, it's important to remember that a
`CliffordNumber{Q,T}` is *not* a concrete type. This is important when creating an `Array` or other
container of `CliffordNumber` elements.

### Internals

A `CliffordNumber{Q,T,L}` is backed by an `NTuple{L,T} where T<:Union{Real,Complex}`. The
coefficients, however, are not indexed in grade order as is done canonically in most resources.

!!! danger Read that again: `CliffordNumber` indexing is not done in grade order.

Instead, the coefficients are arranged in a binary counted fashion, which allows for better SIMD
optimization.

### Constructing a Clifford number

The inner constructor for `CliffordNumber` is `CliffordNumber{Cl,T,L}(x)`, where `x` is any type
that can be converted to an `NTuple{L,T}`. However, in many cases, the type parameters are
redundant, particularly `L`. For this reason, more constructors exist.

In general, one can use a `Vararg` constructor to directly input the values.
```julia-repl
julia> CliffordNumber{APS}(1, 2, 3, 4, 5, 6, 7, 8)
CliffordNumber{APS,Int}(1, 2, 3, 4, 5, 6, 7, 8)
```
Clifford numbers may also be constructed from real numbers, generating a scalar-valued
`CliffordNumber`:
```julia-repl
julia> CliffordNumber{APS}(1)
CliffordNumber{APS,Int}(1, 0, 0, 0, 0, 0, 0, 0)
```
When constructing a `CliffordNumber` from complex numbers, the type parameters become more
important. By default, it is assumed that the element type of a `CliffordNumber` is a `Real`. If
a complex `CliffordNumber` is desired, this must be stated explicitly.
```julia-repl
julia> CliffordNumber{APS}(1 + im)
CliffordNumber{APS,Int}(1, 0, 0, 0, 0, 0, 0, 1)

julia> CliffordNumber{APS,Complex}(1 + im)
CliffordNumber{APS,Complex{Int}}(1 + im, 0, 0, 0, 0, 0, 0, 0)
```

## Quadratic forms

Before getting started with Clifford numbers, it's important to understand how the dimensionality
of the space is stored. Unlike with other data types such as StaticArrays.jl's `SVector`, the
total number of dimensions in the space is not all the information that needs to be stored. Each
basis vector of the space may square to a positive number, negative number, or zero, defining the
quadratic form associated with the Clifford algebra. This information needs to be tracked as a type
parameter for `CliffordNumber`.

To handle this, the `QuadraticForm{P,Q,R}` type is used to store information about the quadratic
form. In this type, `P` represents

!!! note By convention, the `QuadraticForm` type is not instantiated when used as a type parameter
for `CliffordNumber` instances.

CliffordNumbers.jl provides the following aliases for common algebras:

| Algebra    | Alias                  | Note                                           |
| `VGA{D}`   | `QuadraticForm{D,0,0}` | Vanilla/vector geometric algebra               |
| `PGA{D}`   | `QuadraticForm{D,0,1}` | Projective geometric algebra                   |
| `APS`      | `QuadraticForm{3,0,0}` | Algebra of physical space                      |
| `STA`      | `QuadraticForm{1,3,0}` | Spacetime algebra. By default, uses a -+++     |
|            |                        | convention to distinguish it from a conformal  |
|            |                        | geometric algebra.                             |

Currently, an alias for conformal geometric algebras (`CGA{D}`) does not exist, as it requires some
type parameter trickery that hasn't been figured out yet.
