# Performance tips

While CliffordNumbers.jl is intended to be a fast library without much special attention needed to
obtain good performance, there are some ways to ensure you are maximizing performance that may not
be obvious to a first-time user.

## Use literal powers when possible

When exponentiating an `AbstractCliffordNumber`, `x^2` is preferred to `x*x`. This package overloads
`Base.literal_pow` to provide efficient exponentiation in terms of the smallest type needed to
represent the result.

As an example, `KVector{1}`, being a 1-blade, necessarily squares to a scalar, so the return type
of `(x::KVector{1,Q})^2` is always `KVector{0,Q}`. In general, every even power of a 1-blade is a
scalar, and every odd power is a 1-blade. By contrast, `x*x` returns an `EvenCliffordNumber{Q}`. For
small algebras, this distinction may not be significant, but in the 4D case, such as with the
spacetime algebra, `x*x` is 8 times longer than `x^2` if `x` is a `KVector{1}`.

It is also important to note that variable exponents will usually trigger promotion to either
`EvenCliffordNumber` or `CliffordNumber`, depending on the grading of the input. The greatest
benefit comes from exponents that may be evaluated at compile time.

## Know which scalar types work best

In general, `AbstractFloat` scalars are the best choice for the scalar type, followed by `Integer`,
then `Rational`, then `BigInt` or `BigFloat`.

!!! note Fixed point numbers
    Fixed-point numbers haven't been tested with this package, but they are likely as performant as
    their backing type, which is usually an `Integer`.

### Scalars that are not pure bits

The most important examples of scalars that are not bits types are `BigInt` and `BigFloat`. While
these may be necessary to represent certain coefficients accurately, there is a performance penalty
for using them.

When a bits type is used as the scalar, Julia can represent an `AbstractCliffordNumber` as an inline
array of scalars without having to resort to pointers. This provides two major benefits:
  * The type can live on the stack, meaning that no allocations are needed.
  * The type can be stored contiguously in an `Array`.

For types that are not pure bits, the scalars are stored as pointers to data on the heap, which may
or may not be contiguous, and the performance is usually reduced.

### `Rational`

Unlike the machine integer types `Int8`, `Int16`, `Int32`, and `Int64`, or the floating point types
`Float16`, `Float32`, and `Float64`, performing operations with `Rational` types is significantly
slower, since it requires the use of Euclidean division, which is very slow compared to addition and
multiplication.

Specialized kernels for `Rational` products are provided that minimize the number of division
operations needed, but these are still about 100 times slower than performing the same product with
floating point scalars. Improving the speed of these multiplications is a goal for this package, but
do not expect the speed to ever rival that of floating-point multiplication or division.

## Avoid using `CliffordNumber` in most circumstances

When performing geometric operations, you can almost always restrict yourself to `OddCliffordNumber`
or `EvenCliffordNumber` as your largest types. Avoid adding objects with odd and even grades, as
this will trigger promotion to `CliffordNumber`.

If it appears `CliffordNumber` is needed, you may want to ask whether you are using the correct
algebra for your problem. This question is less so about optimization of your code and more so about
how you conceptualize your problem: sometimes it makes more sense to use a larger algebra to get a
better picture of the geometry, and use `EvenCliffordNumber` in the larger algebra as your data
representation.

One circumstance where you may wish to use `CliffordNumber` is when working with idempotents of the
algebra, which may be relevant when working with spinors.

## Know how your CPU's SIMD extensions work

At this point, you've hit hyperoptimization territory. Notes here are intended more so for
developers contributing to this package than for ordinary users, but they are included for the sake
of completeness.

### SIMD width

CliffordNumbers.jl is written to take advantage of CPU SIMD extensions. Usually, these consist of
fixed-size registers of 128, 256, or 512 bits, depending on the implementation.

While Julia can generate efficient code for `Tuple` types with up to 32 elements, it's especially
efficient to use types that fit into one SIMD register. When using an x86 CPU with AVX2, which
provides 256-bit registers, `EvenCliffordNumber{STA,Float32}` is a more performant choice than
`EvenCliffordNumber{STA,Float64}`, since the former consists of 8 32-bit scalars (256 bits total)
and the latter is twice the size (512 bits total).

### SIMD operations

The various products of geometric algebra can be implemented in terms of three common SIMD
operations: add, multiply, and permute. Although these operations are usually available for every
basic real number type, some implementations have gaps which may prevent SIMD extension from being
used.

One notable pitfall is the lack of SIMD multiply instructions for `Int64` on x86 with AVX2 only. It
is possible to perform a 64-bit multiply using AVX2, but these are not available as single
instructions, and Julia's compiler does not generate SIMD code for these operations. On AVX-512 and
AVX10, this instruction is available.
