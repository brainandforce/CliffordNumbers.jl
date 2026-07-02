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
CliffordNumbers.blade_indices_type
CliffordNumbers.CGNBladeIndices
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

## Determinant, trace, and characteristic polynomial

```@docs
CliffordNumbers.det
CliffordNumbers.tr
CliffordNumbers.charpoly_coeffs
```

## Blade and versor structure

```@docs
CliffordNumbers._structure_rtol
CliffordNumbers._euclidean
CliffordNumbers._parity_energy
CliffordNumbers._plucker_ok
CliffordNumbers._probe_vectors
CliffordNumbers._isoclinic_split
CliffordNumbers._bivector_split
```

## Logarithms and square roots of general multivectors

```@docs
CliffordNumbers._study_split
CliffordNumbers._study_log
CliffordNumbers._study_sqrt
CliffordNumbers._even_content
CliffordNumbers._even_magnitude
```

## Return types for operations

```@docs
CliffordNumbers.product_return_type
CliffordNumbers.exponential_type
```
