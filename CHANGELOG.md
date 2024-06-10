# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
  - Missing exports for `regressive` and `∨` (#14 - thanks @ajahraus!)

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

[Unreleased]: https://github.com/brainandforce/CliffordNumbers.jl
[0.1.4]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.4
[0.1.3]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.3
[0.1.2]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.0
