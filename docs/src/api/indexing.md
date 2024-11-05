# Grades and indices

## Grades

```@docs
CliffordNumbers.nonzero_grades
CliffordNumbers.has_grades_of
```

## `BladeIndex`

```@docs
CliffordNumbers.BladeIndex
CliffordNumbers.is_same_blade
```

### Special indices

```@docs
CliffordNumbers.scalar_index
CliffordNumbers.pseudoscalar_index
```

### Tools for implementing mathematical operations

```@docs
Base.reverse(::BladeIndex)
CliffordNumbers.grade_involution(::BladeIndex)
Base.conj(::BladeIndex)
CliffordNumbers.left_complement(::BladeIndex)
CliffordNumbers.right_complement(::BladeIndex)
```

```@docs
CliffordNumbers.grade(::BladeIndex)
CliffordNumbers.sign_of_square
CliffordNumbers.signbit_of_square
CliffordNumbers.nondegenerate_square
CliffordNumbers.sign_of_mult
CliffordNumbers.signbit_of_mult
CliffordNumbers.nondegenerate_mult
Base.:*(::T, ::T) where T<:BladeIndex
CliffordNumbers.has_wedge
```

## `BitIndices` and related types

```@docs
CliffordNumbers.AbstractBitIndices
CliffordNumbers.BitIndices
CliffordNumbers.TransformedBitIndices
```
