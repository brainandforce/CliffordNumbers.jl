# Mathematical operations

## Involutions and duals

```@docs; canonical=false
Base.reverse(::BitIndex)
Base.adjoint(::BitIndex)
CliffordNumbers.grade_involution(::BitIndex)
Base.conj(::BitIndex)
CliffordNumbers.left_complement(::BitIndex)
CliffordNumbers.right_complement(::BitIndex)
```

## Inverses
```@docs
CliffordNumbers.versor_inverse
Base.inv(::AbstractCliffordNumber)
```

## Addition and subtraction

Addition and subtraction integrate seamlessly with those of the Julia Base number types, and no
special documentation is included here.

## Products

### Products with scalars

The standard multiplication and division operations (`*`, `/`, `//`) between Clifford numbers and
scalars behave as expected. `Base.muladd` has been overloaded to take advantage of fma instructions
available on some hardware platforms.

```@docs
Base.muladd(x::BaseNumber, y::T, z::T) where T<:AbstractCliffordNumber
```

### Geometric products

```@docs
Base.:*
CliffordNumbers.:⨼
CliffordNumbers.:⨽
CliffordNumbers.dot
CliffordNumbers.hestenes_dot
CliffordNumbers.:∧
CliffordNumbers.:∨
CliffordNumbers.:×
CliffordNumbers.:⨰
```

### Scalar products and normalization

```@docs
CliffordNumbers.scalar_product
Base.abs2(::AbstractCliffordNumber)
Base.abs(::AbstractCliffordNumber)
CliffordNumbers.normalize
```

## Exponentiation

```@docs
Base.exp(::AbstractCliffordNumber)
CliffordNumbers.exppi
CliffordNumbers.exptau
```
