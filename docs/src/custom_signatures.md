# Custom metric signatures

The pre-defined signature families cover the algebras most code needs:

  * [`VGA`](@ref CliffordNumbers.Metrics.VGA) — vanilla geometric algebras.
  * [`PGA`](@ref CliffordNumbers.Metrics.PGA) — projective geometric algebras.
  * [`CGA`](@ref CliffordNumbers.Metrics.CGA) — conformal geometric algebras.
  * [`LGA`](@ref CliffordNumbers.Metrics.LGA) — Lorentzian geometric algebras.
  * [`Exterior`](@ref CliffordNumbers.Metrics.Exterior) — exterior algebras.

When none of them fit, or when a downstream package wants a signature type that carries its own
meaning, there are two ways to supply a custom metric: a generic
[`Signature`](@ref CliffordNumbers.Metrics.Signature) value, or a new
[`Metrics.AbstractSignature`](@ref CliffordNumbers.Metrics.AbstractSignature) subtype.

## Building a `Signature` value

An orthonormal metric is fully described by which basis 1-blades square to `+1`, `-1`, or `0`, plus
the index of the first basis vector. The generic [`Signature`](@ref CliffordNumbers.Metrics.Signature)
stores exactly this in a few `isbits` fields: a dimension count, a bitmask of the negative-squaring
dimensions, a bitmask of the degenerate (zero-squaring) dimensions, and a first index.

```julia
using CliffordNumbers

# Phase space (2, 2, 0): 2 dimensions squaring to +1, 2 to -1, none degenerate.
# The negative dimensions are 3 and 4, so the negative mask is 0b1100 = 0x0c.
const PhaseSpace2D = Signature(4, 0x0c, 0x00)

x = one(CliffordNumber{PhaseSpace2D,Float64})
```

The value is placed directly in the algebra type parameter (`CliffordNumber{PhaseSpace2D,Float64}`).
Every operation works with it, because the generated kernels read the metric straight from the
value's fields. The package's own projective spacetime algebra
[`STAP`](@ref CliffordNumbers.Metrics.STAP) is defined this way.

Pass a fourth argument to set a non-default first index. Projective and Lorentzian algebras
conventionally start at `0` or below: `Signature(5, 0b11100, 0b00001, -1)` is `STAPWest`.

## Reading the metric with `metric_tuple`

[`metric_tuple`](@ref CliffordNumbers.Metrics.metric_tuple) returns a signature's metric as an
`NTuple{N,Int8}` of the basis 1-blade squares, ordered from `firstindex`:

```jldoctest
julia> using CliffordNumbers

julia> metric_tuple(VGA(3))
(1, 1, 1)

julia> metric_tuple(STA)
(1, -1, -1, -1)

julia> metric_tuple(PGA(2))
(0, 1, 1)
```

This is the single source of metric data inside the package: the bit masks
`positive_square_bits`/`negative_square_bits`/`zero_square_bits` that the multiplication kernels use
are derived from it. It is allocation-free and constant-folded when the signature is a compile-time
constant, which it is whenever the signature lives in a type parameter.

## Defining a signature type

A dedicated type is useful for a family of related algebras or to attach domain meaning to a
signature. Subtype [`Metrics.AbstractSignature`](@ref CliffordNumbers.Metrics.AbstractSignature) and
implement the interface: `dimension`, `firstindex`, and `getindex` are required, while
`is_degenerate` and `is_positive_definite` are recommended.

The phase-space signature `(D, D, R)` — `D` dimensions squaring to `+1`, `D` squaring to `-1`, and
`R` degenerate dimensions — makes a representative example. Both counts are type parameters, so one
type describes the whole family:

```julia
using CliffordNumbers
using CliffordNumbers.Metrics

struct PhaseSpace{D,R} <: Metrics.AbstractSignature end

Metrics.dimension(::PhaseSpace{D,R}) where {D,R} = 2D + R
Base.firstindex(::PhaseSpace) = 1
Base.@propagate_inbounds function Base.getindex(s::PhaseSpace{D}, i::Int) where D
    @boundscheck checkbounds(s, i)
    return i <= D ? Int8(1) : (i <= 2D ? Int8(-1) : Int8(0))
end
Metrics.is_degenerate(::PhaseSpace{D,R}) where {D,R} = R > 0
Metrics.is_positive_definite(::PhaseSpace) = false
```

`dimension` returns the total number of basis 1-vectors, `2D + R`, and `getindex` reports each one's
square; the two agree through the shared `D` and `2D` thresholds. `PhaseSpace{2,0}` is then the
4-dimensional `(2, 2, 0)` algebra and `PhaseSpace{3,1}` the 7-dimensional `(3, 3, 1)` algebra. The
type works with `metric_tuple`, `dimension`, and everything else that reads the metric at run time.

## Generated kernels and world age

CliffordNumbers.jl compiles its hot paths — the geometric product, the automorphisms, the scalar
product, and `charpoly_coeffs` — with
[`@generated` functions](https://docs.julialang.org/en/v1/manual/metaprogramming/#Generated-functions).
A generated function's generator runs in the world age in which the generated function was defined,
which is when CliffordNumbers.jl was loaded. Methods defined afterwards, in a downstream package, are
not visible to it.

A bare subtype value therefore cannot drive these kernels. Inside the generator, dispatch on
`dimension(Q)` and `getindex(Q, i)` falls back to a generic method — the subtype's overrides do not
exist at that world age — and the product throws:

```julia
Q = PhaseSpace{2,0}()
x = CliffordNumber{Q,Float64}(ntuple(_ -> 1.0, 16))
x * x          # throws: the mul generator cannot see PhaseSpace's interface methods
```

Reduce the subtype to a [`Signature`](@ref CliffordNumbers.Metrics.Signature) value with the
one-argument [`Signature(s)`](@ref CliffordNumbers.Metrics.Signature) converter, and parameterize the
Clifford numbers on that value:

```julia
const PS = Signature(PhaseSpace{2,0}())   # → Signature(4, 0x0c, 0x00, 1)

C = CliffordNumber{PS,Float64}
x = C(ntuple(_ -> 1.0, 16))
x * x          # works: the metric is in isbits fields the generator reads directly
```

`Signature(s)` runs at the caller's world age, so it can call the subtype's `dimension` and
`getindex`, and folds them into a value the kernels accept. A constructor that returns the value
directly is a convenient entry point:

```julia
phasespace(d::Integer, r::Integer = 0) = Signature(PhaseSpace{d,r}())
```

Keep the subtype for dispatch, `show`, and domain-specific type parameters, and build the Clifford
numbers over `Signature(s)`.

## Summary

| Goal | Approach |
| --- | --- |
| One-off custom algebra | a [`Signature`](@ref CliffordNumbers.Metrics.Signature) value in the type parameter |
| A typed signature with domain meaning | an `AbstractSignature` subtype, with `Signature(s)` for the Clifford-number parameter |
| Read a signature's metric | [`metric_tuple`](@ref CliffordNumbers.Metrics.metric_tuple) |
