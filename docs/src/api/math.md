# Mathematical operations

## Grade automorphisms

```@docs; canonical=false
Base.reverse(::BitIndex)
CliffordNumbers.grade_involution(::BitIndex)
Base.conj(::BitIndex)
```

## Duals and inverses

```@docs
CliffordNumbers.dual
CliffordNumbers.undual
CliffordNumbers.versor_inverse
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
Base.:*(::AbstractCliffordNumber{Q}, ::AbstractCliffordNumber{Q}) where Q
CliffordNumbers.left_contraction
CliffordNumbers.right_contraction
CliffordNumbers.dot
CliffordNumbers.hestenes_product
CliffordNumbers.wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
CliffordNumbers.commutator
CliffordNumbers.anticommutator
CliffordNumbers.sandwich
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
