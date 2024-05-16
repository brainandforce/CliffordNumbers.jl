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
    In the future, `Base.:~` will no longer be used for this operation; instead `Base.adjoint` will
    be overloaded, providing `'` as a syntax for the reverse.

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

When interpreting even multivectors as elements of the even subalgebra of the algebra of interest,
the grade involution in the even subalgebra is equivalent to the reverse in the algebra of interest.

Grade involution is equivalent to complex conjugation in when dealing with the even subalgebra of 2D
space, which is isomorphic to the complex numbers, but this is *not* true for quaternion
conjugation. Instead, use the Clifford conjugate (described below).

```@docs; canonical=false
CliffordNumbers.grade_involution(::BitIndex)
```

#### Clifford conjugation

The *Clifford conjugate* is the combination of the reverse and grade involution. This is available
via an overload of `Base.conj`.

!!! warning
    `conj(::AbstractCliffordNumber)` implements the Clifford conjugate, not the reverse!

When dealing with the even subalgebras of 2D and 3D VGAs, which are isomorphic to the complex
numbers and quaternions, respectively, the Clifford conjugate is equivalent to complex conjugation
or quaternion conjugation. Otherwise, this is a less widely used operation than the above two.

```@docs; canonical=false
Base.conj(::BitIndex)
```

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
    You can define elements of exterior algebras directly by using `QuadraticForm{0,0,R}`, whose
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
`CliffordNumber.hestenes_product`, which is also not exported.

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

### Defining new products: Multiplication internals

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
