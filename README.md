# CliffordNumbers

[![Stable][docs-stable-img]][docs-stable-url]
[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][ci-status-img]][ci-status-url]
[![Coverage][codecov-img]][codecov-url]

A simple, statically sized multivector (Clifford number) implementation for Julia using graded
representations. This allows for common multivector operations, particularly the various products
of geometric algebra, to be easily implemented with extremely high performance (faster than matrix
multiplications of matrix representations) without depending on any linear algebra library.
Additionally, the multivectors provided by this package can be stored inline in arrays or other data
structures.

# Clifford numbers

## Why use them?

Clifford algebras are algebras of orthonormality: they expand a real or complex vector space with a
notion of normal vectors (those with unit square, or unit magnitude) and orthogonal vectors.

The term "[geometric algebra][ga-wikipedia]" is synonymous in a mathematical sense with Clifford
algebra, but is used by a practicioners of a new pedagogical movement to emphasize their use in
applied mathematics. In particular, *additive representations* are favored over matrix
representations for reasons of clarity, and the package implements all representations of elements
and operations on them in this manner.

If you work with 3D graphics, you may already be familiar with the use of Clifford algebras: the
quaternions used to represent rotation are a Clifford algebra. Specifically, they are isomorphic to
the even subalgebra of the [algebra of physical space][aps-wikipedia], and the elements represent
rotations in 3D. The full algebra of physical space augments rotations with reflections, allowing
any point isometry (combined with dilations) to be represented.

Clifford algebras are also ubiquitous in physics, as the orthonormality relationship applies to 3D
space and (3+1)D spacetime. The Pauli matrices are a matrix representation of the algebra of
physical space, and the Dirac matrices are a matrix representation of the
[spacetime algebra][sta-wikipedia]. These Clifford algebras naturally extend vector algebra, and
even allow for calculus to be done with them.

Here is a short guide to the most commonly used algebras, which can be extended to any number of
dimensions.
  * **VGA (vanilla geometric algebra):** A drop-in replacement for standard vector algebra, suitable
    for classical and quantum physics. Common operations like the cross product and dot product have
    equivalents in VGA, and can be extended to spaces of arbitrary dimension.
  * **PGA (projective geometric algebra):** Extends VGA with a degenerate dimension so that points,
    lines, planes, and related objects can be represented at arbitrary offsets from the origin. Not
    only can it represents point isometries, it can also seamlessly combine them with arbitrary
    translations.
  * **CGA (conformal geometric algebra):** Extends VGA with one positive squaring and one negative
    squaring dimension, and allows for the representation of *k*-spheres. Two-dimensional CGA is
    the algebra of compass and straightedge constructions.
  * **STA (spacetime algebra):** The algebra of Minkowski space, with spatial dimensions squaring to
    values of the opposite sign of temporal dimensions.

Elements of all of the above algebras are constructible in this package.

## Types

This package exports `AbstractCliffordNumber{Q,T}` and its subtypes, which describe the behavior of
multivectors with algebra `Q` and scalar type `T<:Union{Real,Complex}`. This is a subtype of
`Number`, and therefore acts as a scalar for the purpose of broadcasting. For this reason, we
provide an `nblades` function separate from `Base.length` to count the number of blades represented
by a type.

To index an `AbstractCliffordNumber`, we provide the `BitIndex{Q}` type, which allows for arbitrary
components to be indexed, and the `BitIndices{Q,C<:AbstractCliffordNumber{Q}}` type, which provides
all valid indices of instances of `C` that are not constrained to be zero. Indexing with ordinary
integers is disallowed, but `Tuple(::AbstractCliffordNumber)` obtains the backing `Tuple`.

`AbstractCliffordNumber{Q,T}` includes the following concrete subtypes:
  * `CliffordNumber{Q,T,L}`, which represents the coefficients associated with all basis blades.
  * `EvenCliffordNumber{Q,T,L}` and `OddCliffordNumber{Q,T,L}`, which represents multivectors with
only basis blades of even or odd grade being nonzero. These are especially important when dealing
with physically realizable Euclidean transformations (rotations and translations).
  * `KVector{K,Q,T,L}`, which represents multivectors with only basis blades of grade `K` being
nonzero. This is especially useful for representing common vectors, bivectors, pseudovectors, etc.

## Promotion and conversion

The type promotion system is heavily leveraged to minimize the memory footprint needed to represent
the results of various operations. Promotion can convert the numeric types associated with two 
`AbstractCliffordNumber{Q}` instances (for instance, the sum of `KVector{1,APS,Int}` and
`KVector{1,APS,Float64}` is `KVector{1,APS,Float64}`), but it can also leverage grade information to
promote to smaller types: the sum of `KVector{1,APS,Int}` and `KVector{3,APS,Int}` is an
`OddCliffordNumber{APS,Int}`, but the sum of `KVector{1,APS,Int}` and `KVector{2,APS,Int}` is a
`CliffordNumber{APS,Int}`.

We provide the `scalar_promote` and `scalar_convert` functions to allow for promotion of the scalar
types backing an `AbstractCliffordNumber` without needlessly expanding the represented grades.

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
  * Addition (`+`), subtraction and negation (`-`) between algebra elements, or between algebra
    elements and scalars
  * Products of the Clifford algebra:
      * The geometric product  (`*`)
      * The wedge product (`∧`)
      * Left (`⨼`) and right (`⨽`) contractions
      * The dot product and the Hestenes dot product
      * The commutator product (`×`) and anticommutator product (`⨰`)
  * Scalar left (`/`) and right (`\`) division, including rational division (`//`)
  * Efficient `muladd` operations involving scalars and multivectors
  * The reverse (`'`), grade involution, and Clifford conjugation
  * The modulus and absolute value (with `abs2` and `abs`)
  * Exponentiation
      * Efficiently raising multivectors to integer powers
      * Natural exponentials of multivectors

Some of the names or implementations of various operations may change in the near future.

Note that `AbstractCliffordNumber{Q,T}` is a scalar type, so dotted operators do not apply any
operations to each coefficient individually. However, they can be used to perform elementwise
operations on collections of Clifford numbers.

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
[aps-wikipedia]:    https://en.wikipedia.org/wiki/Algebra_of_physical_space
[sta-wikipedia]:    https://en.wikipedia.org/wiki/Spacetime_algebra
[ga-wikipedia]:     https://en.wikipedia.org/wiki/Geometric_algebra
