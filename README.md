# CliffordNumbers

[![Stable][docs-stable-img]][docs-stable-url]
[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][ci-status-img]][ci-status-url]
[![Coverage][codecov-img]][codecov-url]

A simple, fully static multivector implementation for Julia. While not the most space efficient, it
allows for fast prototyping and implementation of geometric algebras and multivectors of arbitrary
dimension or metric signature.

# Clifford numbers

## Types

This package exports `AbstractCliffordNumber{Q,T}` and its subtypes, which describe the behavior of
multivectors with quadratic form `Q` and scalar type `T<:Union{Real,Complex}`. This is a subtype of 
`Number`, and therefore acts as a scalar.

`AbstractCliffordNumber{Q,T}` includes the following concrete subtypes:
  * `CliffordNumber{Q,T,L}`, which represents the coefficients associated with all basis blades.
  * `EvenCliffordNumber{Q,T,L}` and `OddCliffordNumber{Q,T,L}`, which represents multivectors with
only basis blades of even or odd grade being nonzero. These are especially important when dealing
with physically realizable Euclidean transformations (rotations and translations).
  * `KVector{K,Q,T,L}`, which represents multivectors with only basis blades of grade `K` being
zero. This is especially useful for representing common vectors and bivectors.

## Promotion

The type promotion system is heavily leveraged to minimize the memory footprint needed to represent
the results of various operations. Promotion can convert the numeric types associated with two 
`AbstractCliffordNumber{Q}` instances (for instance, the sum of `KVector{1,APS,Int}` and
`KVector{1,APS,Float64}` is `KVector{1,APS,Float64}`), but it can also leverage grade information to
promote to smaller types: the sum of `KVector{1,APS,Int}` and `KVector{3,APS,Int}` is an
`OddCliffordNumber{APS,Int}`, but the sum of `KVector{1,APS,Int}` and `KVector{2,APS,Int}` are
`CliffordNumber{APS,Int}`.

## Indexing

Although `AbstractCliffordNumber` instances are scalars, the `BitIndex{Q}` type can be used to
retrieve coefficients associated with specific basis blades. The full set of `BitIndex{Q}` types for
some `x::AbstractCliffordNumber` can be generated with `BitIndices(x)`, and this is a binary ordered
vector of `BitIndex{Q}` objects.

Mathematical operations are defined generically by working with the `BitIndex{Q}` objects associated
with an `AbstractCliffordNumber{Q,T}`. Elementwise operations on each element of a `BitIndices`
instance returns a `TransformedBitIndices`, a wrapper which lazily associates a function with a
`BitIndices` object, and this can be used to implement grade dependent operations, such as
(anti)automorphisms.

## Operations

The following mathematical operations are supported by this package:
  * Addition (`+`), subtraction and negation (`-`)
  * The geometric product (`*`)
  * Scalar left (`/`) and right (`\`) division, including rational division (`//`)
  * The reverse (`'`), grade involution, and Clifford conjugation
  * The modulus and absolute value (with `abs2` and `abs`)
  * The wedge product (`∧`)
  * The left (`⨼`) and right (`⨽`) contractions
  * The dot product and the Hestenes dot product
  * The commutator product (`×`)
  * Exponentiation

Some of the names or implementations of various operations may change in the near future.

Note that `AbstractCliffordNumber{Q,T}` is a scalar type, so dotted operators do not apply any
operations to each coefficient individually. However, they can be used to perform elementwise
operations on collections of Clifford numbers.

# Features to be added or improved

## Type system

  * A `StaticMultivector{Q,T,L}` type, allowing for non-static implementations.
  * `numeric_type` will likely be renamed to `scalartype`.

## Relationships between Clifford algebras

  * Determining the even subalgebra associated with some Clifford algebra.
  * Converting even multivectors of an algebra to multivectors in the even subalgebra.

## Mathematical operations

  * Better numerical accuracy and stability for some operations. In particular, properly leveraging
`sinpi`, `cospi`, `hypot`, and other methods with better numerical accuracy and stability 
internally.

## Interoperability

  * How do we handle the dot product of this package and the dot product of the `LinearAlgebra`
standard library? This package has no dependencies, but the semantics should be compatible as much
as possible.
  * Secondarily, `StaticArrays.similar_type` and `CliffordNumbers.similar_type` have essentially
identical behavior.

[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://brainandforce.github.io/CliffordNumbers.jl/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://brainandforce.github.io/CliffordNumbers.jl/dev
[ci-status-img]:    https://github.com/brainandforce/CliffordNumbers.jl/workflows/CI/badge.svg
[ci-status-url]:    https://github.com/brainandforce/CliffordNumbers.jl/actions
[aqua-img]:         https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]:         https://github.com/JuliaTesting/Aqua.jl
[codecov-img]:      https://codecov.io/gh/brainandforce/CliffordNumbers.jl/branch/main/graph/badge.svg
[codecov-url]:      https://codecov.io/gh/brainandforce/CliffordNumbers.jl/
