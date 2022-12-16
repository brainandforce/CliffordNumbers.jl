# CliffordNumbers

[![Stable][docs-stable-img]][docs-stable-url]
[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][ci-status-img]][ci-status-url]
[![Coverage][codecov-img]][codecov-url]

A simple, fully static multivector implementation for Julia. While not the most space efficient, it
allows for fast prototyping and implementation of geometric algebras and multivectors of arbitrary
dimension or metric signature.

Currently, the package exports the `CliffordNumber{Cl,T}` type, which represents an element of the
Clifford algebra `Cl` with elements of type `T`, which have `L` elements. Future releases will
include more storage-efficient data structures for common, sparser elements such as k-vectors and
blades.

Currently supported operations are addition (`+`), subtraction and negation (`-`), the geometric
product (`*`), and the Hodge star (`⋆`).

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
