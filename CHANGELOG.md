# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
[0.1.2]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/brainandforce/CliffordNumbers.jl/releases/tag/v0.1.0
