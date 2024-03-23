# Internal methods

## Hamming weights and related operations

```@docs
CliffordNumbers.isevil
CliffordNumbers.isodious
CliffordNumbers.number_of_parity
CliffordNumbers.evil_number
CliffordNumbers.odious_number
CliffordNumbers.next_of_hamming_weight
CliffordNumbers.hamming_number
```

## Indexing

```@docs
CliffordNumbers.signmask
CliffordNumbers._sort_with_parity!
```

## Construction

```@docs
CliffordNumbers.zero_tuple
CliffordNumbers.check_element_count
```

## Multiplication kernels

```@docs
CliffordNumbers.mul
CliffordNumbers.GradeFilter
CliffordNumbers.nondegenerate_mask
CliffordNumbers.mul_mask
CliffordNumbers.bitindex_shuffle
CliffordNumbers.widen_grade_for_mul
```

## Taylor series exponential
```@docs
CliffordNumbers.exp_taylor
```

## Return types for operations

```@docs
CliffordNumbers.product_return_type
CliffordNumbers.exponential_type
```

## Unused features

```@docs
CliffordNumbers.RepresentedGrades
```
