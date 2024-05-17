```@meta
CurrentModule = CliffordNumbers
```

# CliffordNumbers

[CliffordNumbers.jl](https://github.com/brainandforce/CliffordNumbers.jl) is a package that
provides fully static multivectors (Clifford numbers) in arbitrary dimensions and metrics. While
in many cases, sparse representations of multivectors are more efficient, for spaces of low
dimension, dense static representations may provide a performance and convenience advantage.

## Design goals

The goal of this package is to provide a multivector implementation that:
- Allows for the construction of multivectors in arbitrary metrics, with coefficients that subtype any instance of Julia's base numeric types, `Real` or `Complex`.
- Provides data structures of fixed sizes that represent multivectors. This allows for instances to be allocated on the stack or stored inline in an `Array` rather than as pointers to individually allocated instances.
- Provides dense representations of multivectors, as well as convenient sparse representations, which can be constructed from each other, converted in a way that guarantees representability, and allows for promotion between instances.
- Subtypes `Number`: The term "Clifford number" emphasizes the perspective of multivectors as an extension of the real numbers, in the same way that complex numbers and quaternions extend them. (It should be noted that both complex numbers and quaternions are Clifford algebras themselves!)
- Aggressively optimizes all mathematical operations, utilizing `fma` operations and SIMD instructions whenever possible.
- Interoperates with automatic differentiation tools and other packages which allow for the implementation of operations from geometric calculus.
