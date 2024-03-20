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
