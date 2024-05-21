---
title: 'CliffordNumbers.jl: '
tags:
  - Julia
  - geometric algebra
  - Clifford algebra
authors:
  - name: Brandon S. Flores
    orcid: 0000-0003-0816-4491
    affiliation: 1
    corresponding: true
affiliations:
 - name: Department of Chemistry, Univesity of Wisconsin-Madison, USA
   index: 1
date: 19 May 2024
bibliography: paper.bib

---

# Summary

For decades, linear algebra libraries have been developed and intensely optimized to maximize
performance for a variety of important operations, including matrix multiplication and numerous
matrix factorizations. Although geometric algebra has gained popularity, particularly in the last
20 years, its development on the software side has not been driven by an overarching shared need in
the way linear algebra libraries have, instead being pushed forward by significantly smaller groups
of both researchers and enthusiasts.

`CliffordNumbers.jl` provides a high-performance implementation of commonly used operations in
geometric algebras of up to 6 dimensions without relying on a linear algebra library. In common
cases, the performance of this library allows for implementations of isometry composition or 
application with significantly greater speed

# Statement of need

Clifford algebras are ubiquitous throughout a number of domains, in particular the physical sciences
and computer graphics. In essence, they are algebras of orthonormality, and can be used to model
isometries: rotations, reflections, translations, and more. Although Clifford algebras are
well-studied mathematical objects, they are usually approached as matrix representations; most
notably the Pauli matrices of quantum mechanics and the gamma matrices found in the Dirac equation.

In recent years, working with Clifford algebras by using an *additive represenatation* not requiring
matrices has gained popularity, especially in practical applications. Although this approach was
initially developed by William Kingdon Clifford himself, it lost out to the vector algebra formalism
of Josiah Willard Gibbs and Oliver Heaviside, and was nearly forgotten until it was revisited by
David Hestenes, who described them as *geometric algebras*. In this article, the terms *Clifford 
algebra* and *geometic algebra* refer to identical mathematical structures, but we will use the term
*geometric algebra* to refer to the study of these structures to model geometry for practical
applications.

Even without explicitly acknowledging geometric algebras, we still rely on them: the Pauli spin
matrices of quantum mechanics are a matrix representation of the *algebra of physical space* (APS),
the Clifford algebra $\text{Cl}_{3,0,0}\left(\mathbb{R}\right)$. In computer graphics, quaternions,
which are used to compose and apply rotations, are isomorphic to the *rotors* of APS, or elements
which represent rotations. Like with quaternions, single-sided multiplications constitute 
composition of rotations, and double-sided multiplications represent application of the rotation to
some object of the space, and this is always faster than an equivalent matrix-vector multiplication
performing the rotation.

# Implementation

`CliffordNumbers.jl` provides a generic implementation of geometric algebras of dimension only
limited by the machine it runs on. The `AbstractCliffordNumber{Q,T}` abstract type is a subtype of
the Julia Base `Number` type, and much of the semantics of the concrete subtypes of
`AbstractCliffordNumber{Q,T}` match those of other `Number` subtypes, including `Real`, 
`Complex{T}`, and the `Quaternion{T}` type provided by `Quaternions.jl`. This design decision was
made because the complex numbers and quaternions are both isomorphic to real Clifford algebras, and

This library does not depend on a linear algebra package. Unlike many other software implementations
of geometric algebras, this package does not use matrix representations to calculate the geometric
product or any related products. Instead, it uses the additive representation internally, where the
coefficients of an `AbstractCliffordNumber{Q,T}` instance are directly associated with the basis
blades of the algebra. Suprisingly, this allows for significant performance improvements over other
implementations that use matrix representations.

# Indexing

Although `AbstractCliffordNumber{Q,T}` instances behave like scalars, it is often necessary to
obtain individual components of an instance, either in the implementation of 

# Performance

For a worst-case comparison, we will provide a benchmark that compares the geometric product of two
arbitrary multivectors of APS, which have 8 elements each, with the usual implementation as a
multiplication of two 2×2 complex-valued matrices. The code for each result is provided, and these
comparisons were performed on an Intel Core i5-13600K with 32 GB RAM:

```
julia> @benchmark @inline $m1 * $m2
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  2.548 ns … 208.505 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.561 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.596 ns ±   2.104 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

              ▂ ▃ ▅ ▆ █ ▇▇ █ █ ▆ ▅ ▃ ▂                         
  ▂▂▁▂▁▃▁▄▁▅▁▆█▁█▁█▁█▁█▁██▁█▁█▁█▁█▁█▁██▁▆▁▅▁▄▁▃▁▃▂▁▂▁▂▁▂▁▂▁▂▂ ▄
  2.55 ns         Histogram: frequency by time        2.58 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark @inline $c1 * $c2
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  2.341 ns … 25.916 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.358 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.378 ns ±  0.434 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

       ▇▅█▂▁                                                  
  ▂▂▃▆██████▄▄▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▁▁▂▂▂▁▂▂▁▁▁▁▁▂▂▂▂▁▂▂▂ ▃
  2.34 ns        Histogram: frequency by time        2.49 ns <

 Memory estimate: 0 bytes, allocs estimate: 0
```

To understand this result, we can analyze the algorithm for matrix multiplication and compare it to
the implementation of the geometric product provided by this package. We want to compare the total
number of operations on real numbers, so we must break down complex addition and multiplication into
its constituent real additions and multiplications. Complex addition is the straightforward addition
of real and imaginary coefficients, so it is two real additions:
```julia
+(z::Complex, w::Complex) = Complex(real(z) + real(w), imag(z) + imag(w))
```
The multiplication of two complex numbers consists of 4 multiplications and 2 additions:
```julia
*(z::Complex, w::Complex) = Complex(real(z) * real(w) - imag(z) * imag(w),
                                    real(z) * imag(w) + imag(z) * real(w))
```
The multiplication of a pair of 2×2 matrices consists of 4 sets of dot products between rows of the
first matrix and columns of the second. These dot products consist of two complex multiplications
and a complex addition, so in total there are 8 complex multiplications and 4 complex additions,
which total 32 real multiplications and 24 real additions.

The geometric product in the additive representation can be implemented with a multiplication table,
requiring 56 multiplications and 56 additions.

However, in practice, isometries of the space are represented by $k$-versors, the geometric product
of $k$ 1-blades. These necessarily only contain elements of either even or odd grades, meaning that
only half as many coefficients are present. On top of this, the `KVector{K,Q,T,L}` data type can be
used to represent $k$-blades, and these only have `binomial(dimension(Q), K)` elements

# Examples

# Development roadmap

The next major development milestones for `CliffordNumbers.jl` will include 

# Dedication and acknowledgements

The author dedicates this paper to Dr. Michael E. Davies and Kristel M. Forlano.

The author acknowledges the community at [bivector.net](https://www.bivector.net), especially 
members of the affiliated Discord server, in contributing their collective knowledge of geometric
algebras, their applications, and software implementations. In particular he would like to
specifically acknowledge the following community members:
- sudgylacmoe, YouTuber and creator of the video A Swift Introduction to Geometric Algebra.
- Joseph Wilson (Jollywatt), author of another geometric algebra package.

# References
