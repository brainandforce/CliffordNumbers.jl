# Clifford numbers

## Supertype and associated functions

```@docs
CliffordNumbers.AbstractCliffordNumber
CliffordNumbers.nblades
CliffordNumbers.scalar_type
CliffordNumbers.similar_type
```

## Concrete types

```@docs
CliffordNumbers.CliffordNumber
CliffordNumbers.Z2CliffordNumber
CliffordNumbers.EvenCliffordNumber
CliffordNumbers.OddCliffordNumber
CliffordNumbers.KVector
CliffordNumbers.grade(::Type{<:KVector{K}}) where K
```

## Promotion and conversion

```@docs
CliffordNumbers.scalar_convert
CliffordNumbers.scalar_promote
Base.widen(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T}
CliffordNumbers.widen_grade
```

## Real and complex algebras

```@docs
Base.real(::AbstractCliffordNumber)
Base.complex(::AbstractCliffordNumber)
```

## Scalar and pseudoscalar components

```@docs
CliffordNumbers.isscalar
CliffordNumbers.ispseudoscalar
CliffordNumbers.scalar
CliffordNumbers.pseudoscalar
```
