# Operations

Like with other numbers, standard mathematical operations are supported that relate Clifford numbers
to elements of their scalar field and to each other.

## Unary operations

### Grade automorphisms

Grade automorphisms are operations which preserves the grades of each basis blade, but changes their
sign depending on the grade. All of these operations are their own inverse.

All grade automorphisms are applicable to `BitIndex` objects, and the way they are implemented is
through constructors that use `TransformedBitIndices` objects to alter each grade.

#### Reverse

The *reverse* is an operation which reverses the order of the wedge product that constructed each
basis blade. This is implemented with methods for `Base.reverse` and `Base.:~`.

!!! note "Syntax changes"
    `Base.:~` for the reverse is deprecated and will be removed. `Base.adjoint` provides `'` as a
    syntax for the reverse, and is the preferred method for performing the reverse.

This is the most commonly used automorphism, and in a sense can be thought of as equivalent to
complex conjugation. When working with even elements of the algebras of 2D or 3D space, this
behaves identically to complex conjugation and quaternion conjugation. However, this is *not* the
case when working in the even subalgebras.

```@docs; canonical=false
Base.reverse(::BitIndex)
```

#### Grade involution

*Grade involution* changes the sign of all odd grades, an operation equivalent to mirroring every
basis vector of the space. This can be acheived with the `grade_involution` function.

When interpreting even multivectors as elements of the even subalgebra of a given algebra, the
reverse operation of the algebra is equivalent to grade involution in the even subalgebra.

Grade involution is equivalent to complex conjugation in when dealing with the even subalgebra of 2D
space (which is isomorphic to the complex numbers), but this is *not* true for quaternion
conjugation. Instead, use the Clifford conjugate (described below).

```@docs; canonical=false
CliffordNumbers.grade_involution(::BitIndex)
```

#### Clifford conjugation

The *Clifford conjugate* is the combination of the reverse and grade involution. This is implemented
as `Base.conj(::AbstractCliffordNumber)`. This operation arises in the application of a
transformation, because the grade involution accounts for the sign change associated with the parity
of the isometry, and it is combined with the reverse to perform the final operation.

!!! warning
    `conj(::AbstractCliffordNumber)` implements the Clifford conjugate, not the reverse!

When dealing with the even subalgebras of 2D and 3D VGAs, which are isomorphic to the complex
numbers and quaternions, respectively, the Clifford conjugate is equivalent to complex conjugation
or quaternion conjugation. Otherwise, this is a less widely used operation than the above two.

```@docs; canonical=false
Base.conj(::BitIndex)
```

#### Inverse

Many elements of a Clifford algebra have an inverse. In general, the inverse of a rotor `R` is equal
to `R' / abs2(R)`. It is possible to use other methods to find inverses of arbitrary multivectors of
certain Clifford algebras, but this has not been implemented yet.

Inverses *cannot* exist for multivectors which square to zero. This is trivially true in any 
Clifford algebra with a degenerate metric, but it is possible to encounter these kinds of elements
in algebras with positive-definite metrics.

## Binary operations

### Addition and subtraction

Addition and subtraction work as expected for Clifford numbers just as they do for other numbers.
The promotion system handles all cases where objects of mixed type are added.

### Products

Clifford algebras admit a variety of products. Common ones are implemented with infix operators.

#### Geometric product

The geometric product, or Clifford product, is the defining product of the Clifford algebra. This is
implemented with the usual multiplication operator `*`, but it is also possible to use parenthetical
notation as it is with real numbers.

#### Wedge product

The wedge product is the defining product of the exterior algebra. This is available with the
`wedge()` function, or with the `∧` infix operator.

!!! tip
    You can define elements of exterior algebras directly by using `Metrics.Exterior(D)`, whose
    geometric product is equivalent to the wedge product.

#### Contractions and dot products

The contraction operations generalize the dot product of vectors to Clifford numbers. While it is
possible to define a symmetric dot product (and one is provided in this package), the generalization
of the dot product to Clifford numbers is naturally asymmetric in cases where the grade of one
input blade is not equal to that of the other.

For Clifford numbers `x` and `y`, the *left contraction* `x ⨼ y` describes the result of projecting
`x` onto the space spanned by `y`. If `x` and `y` are homogeneous in grade, this product is equal to
the geometric product if `grade(y) ≥ grade(x)`, and zero otherwise. For general multivectors, the
left contraction can be calculated by applying this rule to the products of their basis blades.

The analogous right contraction is only nonzero if `grade(x) ≥ grade(y)`, and it can be calculated
with `⨽`.

The *dot product* is a symmetric variation of the left and right contractions, and provides a looser
constraint on the basis blades: `grade(CliffordNumbers.dot(x,y))` must equal
`abs(grade(x) - grade(y))`. The  *Hestenes dot product* is equivalent to the dot product above, but
is zero if either `x` or `y` is a scalar.

!!! note
    Currently, the dot product is implemented with the unexported function `CliffordNumbers.dot`.
    This package does not depend on LinearAlgebra, so there would be a name conflict if this method
    were exported and both this package and LinearAlgebra were loaded.

Contractions are generally favored over the dot products due to their nicer implementations and
properties, which have fewer exceptions. It is generally recommended that the Hestenes dot product 
be avoided, though it is included in this library for the sake of completeness as
`CliffordNumber.hestenes_dot`, which is also not exported.

#### Commutator and anticommutator products

The *commutator product* (or *antisymmetric product*) of Clifford numbers `x` and `y`, denoted
`x × y`, is equal to `1//2 * (x*y - y*x)`. This product is nonzero if the geometric product of `x` 
and `y` does not commute, and the value represents the degree to which they fail to commute.

The commutator product is the building block of Lie algebras; in particular, the commutator products
of bivectors, which are also bivectors. With the bivectors of 3D space, the Lie algebra is
equivalent to that generated by the cross product, hence the `×` notation.

The analogous *anticommutator product* (or *symmetric product*) is `1//2 * (x*y + y*x)`. This uses
the `⨰` operator, which is not an operator generally used for this purpose, but was selected as it
looks similar to the commutator product, with the dot indicating the similarity with the dot
product, which is also symmetric.

#### Defining new products: Multiplication internals

Products are implemented with the fast multiplication kernel `CliffordNumbers.mul`, which accepts
two Clifford numbers with the same scalar type and a `CliffordNumbers.GradeFilter` object. This
`GradeFilter` object defines a method that takes two or more `BitIndex` objects and returns `false`
if their product is constrained to be zero.

`CliffordNumbers.mul` requires that the coefficient types of the numbers being multiplied are the
same. Methods which leverage `CliffordNumbers.mul` should promote the coefficient types of the
arguments to a common type using `scalar_promote` before passing them to the kernel. Any further
promotion needed to return the final result is handled by the kernel.

In general, it is also strongly recommended to promote the types of the arguments to
`CliffordNumbers.Z2CliffordNumber` or `CliffordNumber` for higher performance. Currently, the
implementation of `CliffordNumbers.mul` is asymmetric, and does not consider which input is longer.
Even in the preferred order, we find that `KVector` incurs a significant performance penalty.

### Exponentiation

Exponentiation can be done using either a Clifford number as a base and an integer exponent,
corresponding to iterations of the geometric product, or a Clifford number may be the exponent
associated with a scalar.

#### Integer powers of Clifford numbers

As with all other `Number` instances, this can be done with the `^` infix operator:
```
julia> k = KVector{1,VGA(3)}(4,2,0)
3-element KVector{1, VGA(3), Int64}:
4e₁ + 2e₂

julia> k^2
1-element KVector{0, VGA(3), Int64}:
20
```
However, it may be worth noting that the types may seem to vary based on what Clifford number is
being exponentiated. The bivector `l` squares to an `EvenCliffordNumber` rather than a `KVector{0}`:
```
julia> l = KVector{2,VGA(3)}(0, 6, 9)
3-element KVector{2, VGA(3), Int64}:
6e₁e₃ + 9e₂e₃

julia> l^2
4-element EvenCliffordNumber{VGA(3), Int64}:
-117
```
This is to be expected: only ``k``-blades (``k``-vectors that are the wedge product of ``k``
1-vectors) can be guaranteed to square to a scalar. In dimensions less than 4, all ``k``-vectors are
``k``-blades, but this cannot be assumed in general: the classic example is ``e_1 e_2 + e_3 e_4`` in
4D VGA, whose square contains a 4-vector term. However, in all dimensions, 0-vectors (scalars) and
1-vectors are guaranteed to be blades, and therefore square to scalars. This is also true of 
pseudoscalars and pseudovectors (``n``-blades and ``(n-1)``-blades in an ``n``-dimensional algebra),
but these cases are not recognized yet.

Furthermore, the types will differ if we exponentiate using a variable rather than a literal 
integer:
```
julia> z = 2
2

julia> k^z
8-element CliffordNumber{VGA(3), Float64}:
20.0

julia> l^z
4-element EvenCliffordNumber{VGA(3), Float64}:
-117.0
```
Although the results compare as equal with `==`, they do not with `===`. The reason for this is type
stability, and the mechanisms Julia provides to reduce unnecessary type conversion when enough
information is known at compile time.

Exponentiation is an inherently problematic operation with regards to type stability. This is simply
illustrated when raising integers to integer powers: if the exponent is a positive integer, or zero,
the result should be an integer, but if the exponent is negative, the result cannot be expressed as
an integer unless the number being exponentiated is a unit.

It may seem like we have to promote every Clifford number that isn't a `KVector{0}` to either
`EvenCliffordNumber` or `CliffordNumber`, but Julia gives us a way to work around this. When any
exponentation occurs with a literal number, Julia replaces the expression with `Base.literal_pow`,
and the exponentiation can be converted to a series of multiplications. This package provides these
definitions so that all Clifford numbers can be efficiently exponentiated. In particular, we deal
with the cases of `KVector{0}` and `KVector{1}` explicitly, so that all of their even powers result
in a `KVector{0}`, and the odd powers of a `KVector{1}` are also a `KVector{1}`. In the future, we
will extend these optimizations to pseudoscalars and pseudovectors.

#### Natural exponentiation

It is also possible to raise a real number to a Clifford number power with the `exp` function, or by
raising the irrational constant `ℯ` to an exponent with `^`.
```
julia> exp(KVector{2,VGA(3)}(pi/2, 0, 0))
4-element EvenCliffordNumber{VGA(3), Float64}:
6.123233995736766e-17 + 1.0e₁e₂

julia> e^KVector{2,VGA(3)}(pi/2, 0, 0)
4-element EvenCliffordNumber{VGA(3), Float64}:
6.123233995736766e-17 + 1.0e₁e₂
```
The most common use case for this operation is to construct rotors from bivectors, as illustrated
above.

Internally, exponentiation is done by a Taylor expansion in the general case, but it is also
possible to simplify the exponetiation of `KVector` instances by identifying the sign of the square
of the input. 

Julia provides the `sinpi`, `cospi`, and related functions that allow for the calculation of 
`sin(pi*x)` or `cos(pi * x)`, respectively, with greater accuracy, especially for large `x`.
Although there is no `exppi` function provided by Julia Base, this packages provides one for use
with Clifford numbers, and the accuracy of exponentiation can be expected to be better for 
negative-squaring blades.
```
julia> exppi(KVector{2,VGA(3)}(1/2, 0, 0))
4-element EvenCliffordNumber{VGA(3), Float64}:
1.0e₁e₂

julia> exptau(KVector{2,VGA(3)}(1/4, 0, 0))
4-element EvenCliffordNumber{VGA(3), Float64}:
1.0e₁e₂
```
This package also provides `exptau`, which calculates `exp(2pi * x)`.
