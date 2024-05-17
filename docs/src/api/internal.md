# Internal methods

## Hamming weights and related operations

```@docs
CliffordNumbers.Hamming
CliffordNumbers.Hamming.isevil
CliffordNumbers.Hamming.isodious
CliffordNumbers.Hamming.number_of_parity
CliffordNumbers.Hamming.evil_number
CliffordNumbers.Hamming.odious_number
CliffordNumbers.Hamming.next_of_hamming_weight
CliffordNumbers.Hamming.hamming_number
```

## Indexing

```@docs
CliffordNumbers.signmask
CliffordNumbers._sort_with_parity!
CliffordNumbers.bitindices_type
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
CliffordNumbers.mul_signs
CliffordNumbers.bitindex_shuffle
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
