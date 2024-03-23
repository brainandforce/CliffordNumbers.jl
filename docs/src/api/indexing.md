# Grades and indices

## Grades

```@docs
CliffordNumbers.nonzero_grades
CliffordNumbers.has_grades_of
```

## `BitIndex`

```@docs
CliffordNumbers.BitIndex
CliffordNumbers.is_same_blade
```

### Special indices

```@docs
CliffordNumbers.scalar_index
CliffordNumbers.pseudoscalar_index
```

### Tools for implementing mathematical operations

```@docs
Base.reverse(::BitIndex)
CliffordNumbers.grade_involution(::BitIndex)
Base.conj(::BitIndex)
```

```@docs
CliffordNumbers.grade(::BitIndex)
CliffordNumbers.sign_of_square
CliffordNumbers.signbit_of_square
CliffordNumbers.nondegenerate_square
CliffordNumbers.sign_of_mult
CliffordNumbers.signbit_of_mult
CliffordNumbers.nondegenerate_mult
Base.:*(::T, ::T) where T<:BitIndex
CliffordNumbers.has_wedge
```

## `BitIndices` and related types

```@docs
CliffordNumbers.AbstractBitIndices
CliffordNumbers.BitIndices
CliffordNumbers.TransformedBitIndices
```
