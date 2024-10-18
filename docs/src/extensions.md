# Extensions

CliffordNumbers.jl provides some extensions to allow for interoperability with other packages which
we anticipate to be commonly used alongside it.

!!! note
    Extensions require Julia 1.9; support for previous versions will be dropped with the 0.2.0
    release of CliffordNumbers.jl.

## [Unitful.jl]

[Unitful.jl] provides support for quantities with associated units.

To make a Clifford number into a unitful quantity, simply multiply it by a unit or other unitful
quantity:
```julia-repl
julia> KVector{1, VGA(3)}(1, 2, 3)u"m"
(1e₁ + 2e₂ + 3e₃) m

julia> KVector{1, VGA(3)}(1, 2, 3) * 1.5u"m"
(1.5e₁ + 3.0e₂ + 4.5e₃) m
```

!!! warn "Constructing Clifford numbers from `Quantity` coefficients fails"
    Attempts to construct a Clifford number from `Quantity` objects will fail, because `Quantity`
    does not subtype `Real` or `Complex`:
    ```julia-repl
    julia> KVector{1,VGA(3)}(1u"m", 2u"m", 3u"m")
    ERROR: ArgumentError:
    ...
    ```
    For this reason, Clifford numbers with mixed units cannot be constructed, and this is 
    intentional: `AbstractCliffordNumber` assumes an orthonormal basis.

Supported operations include the geometric product and wedge product.

[Unitful.jl]: https://github.com/PainterQubits/Unitful.jl
