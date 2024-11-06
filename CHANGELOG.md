# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Next]

### Removed
  - **[BREAKING]** `CliffordNumbers.normalize` is no longer exported to avoid a name conflict with
    `LinearAlgebra.normalize`.

## [0.1.10] - 2024-11-06

### Added
  - [Contribution guidelines.](CONTRIBUTING.md)

### Fixed
  - `isinf`, `isnan`, `isreal`, `isinteger`, `iseven`, and `isodd` did not have working methods for
    `AbstractCliffordNumber` arguments.

## [0.1.9] - 2024-10-28

### Added
  - A package extension for interoperability with [Unitful.jl] (only supported on julia 1.9 and up).
    It provides support for wedge products for `Quantity`, including `Quantity{<:Real}` and
    `Quantity{<:AbstractCliffordNumber}` (geometric products were already supported).
  - `nonzero_grades(::Complex)` is defined and returns `0:0`, like `nonzero_grades(::Real)`.
  - `print(::IO, ::AbstractCliffordNumber)` shows a prettier (but not parseable) representation.

### Changed
  - Indexing and equality checking of `BitIndices` is slightly more efficient.

## [0.1.8] - 2024-09-23

### Added
  - Implementations of `flipsign` and `copysign` for `BitIndex`.
  - `+(::BitIndex)` returns the input.
  - Faster multiplication methods for `Rational` coefficients.

### Changed
  - Geometric products of `BitIndex` objects are faster.
  - `scalar_product` is implemented with a generated function independent of `CliffordNumbers.mul`.
  - `normalize` returns `x` if `abs2(x)` is zero.
  - Indexing of Clifford numbers with `Rational` scalar types is significantly faster due to the use
    of `flipsign` instead of multiplication.
  - `KVector` indexing is now faster.

### Fixed
  - `scalar_product` had no method to handle mismatched scalar types.
  - `abs` now returns a `Real` result: it is equal to `sqrt(abs(abs2(x)))`.
  - Indexing `KVector` does not throw an error if the first entry is `1//0`.

## [0.1.7] - 2024-09-11

### Added
  - `Base.literal_pow` definitions for `KVector{1}` that produce `KVector` results.

### Changed
  - Natural exponential are only optimized for `KVector` types whose `BitIndices` uniformly square
    to the same value.

### Fixed
  - `OddCliffordNumber` constructors now fail to operate on scalars unless there is exactly one odd
    element in the algebra.
  - Broken bullets in `AbstractCliffordNumber` docstring.
  - `StackOverflowError` in `KVector{0}` exponentiation.
  - Incorrect natural exponentials in algebras with non-uniform sign signatures.
  - `isequal(x::AbstractCliffordNumber, y::AbstractCliffordNumber)` gives the correct result, no
    longer equal to `x == y` for signed zeros and NaNs.

## [0.1.6] - 2024-06-27

### Added
  - Definition of `signature(::BitIndex{Q})`.
  - Unexported aliases `CliffordNumbers.CliffordScalar{Q,T} === KVector{0,Q,T,1}`, 
    `CliffordNumbers.CliffordVector{Q,T,L} === KVector{1,Q,T,L}`, and
    `CliffordNumbers.CliffordBivector{Q,T,L} === KVector{2,Q,T,L}`.
  - `@basis_vars` macro that defines basis 1-blade variables for an algebra.

### Changed
  - Split up `src/math.jl` to separate files in `src/math/`.

## [0.1.5] - 2024-06-10

### Fixed
  - Missing exports for `regressive` and `∨` (#14 - thanks @ajahraus!)
  - Actually changed the implementations of `scalar_product`, `abs2`, and `abs`.

## [0.1.4] - 2024-06-06

### Added
  - Linear extension of left and right complements to all Clifford numbers.
  - Regressive product (`∨`) definition.

### Changed
  - `scalar_product`, `abs2`, and `abs` are much faster: they just extract the scalar portion of
    the result of `CliffordNumbers.mul(x, y, CliffordNumbers.GradeFilter{:*}())`

### Fixed
  - Incorrect products for elements of algebras with negative-squaring elements.
  - Incorrect complements for elements of algebras with non-positive-definite signatures.
  - Incorrect normalization of division of scalars by Clifford numbers.

## [0.1.3] - 2024-06-04

### Added
  - Left complement (`left_complement()`) and right complement (`right_complement()`) of `BitIndex`
    objects.

### Fixed
  - Incorrect definition of `abs2` for non-positive-definite and degenerate metrics.

## [0.1.2] - 2024-05-31

### Added
  - `Base.float` and `Base.big` definitions for `AbstractCliffordNumber` types and instances.
  - `Base.literal_pow` definitions for `AbstractCliffordNumber` instances raised to constant powers,
    allowing the powers of `KVector`, `EvenCliffordNumber`, or `OddCliffordNumber` to be inferred
    as either `EvenCliffordNumber` or `OddCliffordNumber` depending on the exponent.

### Changed
  - Geometric products involving pseudoscalars (`KVector{K,Q}` where `K === dimension(Q)`) now
    promote to smaller types if possible.
  - `CliffordNumbers.mul` iterates through the indices of the smaller argument, which drastically
    reduces the performance discrepancy when the multiplication arguments are reversed if SIMD
    vectorization is utilized.

## [0.1.1] - 2024-05-28

### Changed
  - `@inline` annotations have been provided for all products.

## [0.1.0] - 2024-05-17

Initial release of CliffordNumbers.jl

[Next]:       https://github.com/brainandforce/CliffordNumbers.jl/tree/next
[Unreleased]: https://github.com/brainandforce/CliffordNumbers.jl
[0.1.10]:     https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.10
[0.1.9]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.9
[0.1.8]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.8
[0.1.7]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.7
[0.1.6]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.6
[0.1.5]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.5
[0.1.4]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.4
[0.1.3]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.3
[0.1.2]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.2
[0.1.1]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.1
[0.1.0]:      https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.0
[Unitful.jl]: https://github.com/PainterQubits/Unitful.jl
