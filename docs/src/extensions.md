# Extensions

CliffordNumbers.jl provides some extensions to allow for interoperability with other packages which
we anticipate to be commonly used alongside it.

## [Quaternions.jl]

[Quaternions.jl] provides the `Quaternion` type.

This extension treats `Quaternion{T}` as identical to `EvenCliffordNumber{VGA(3),T,4}`, and each 
type can be constructed from the other. This also extends to `AbstractCliffordNumber{VGA(3)}`.

```julia-repl
julia> q = Quaternion(0, 1, 2, 3)
Quaternion{Int64}(0, 1, 2, 3)

julia> e = EvenCliffordNumber{VGA(3)}(q)
4-element EvenCliffordNumber{VGA(3), Int64}:
1e₁e₂ + 2e₁e₃ + 3e₂e₃

julia> Quaternion(e)
Quaternion{Int64}(0, 1, 2, 3)

julia> KVector{2,VGA(3)}(q)
3-element KVector{2, VGA(3), Int64}:
1e₁e₂ + 2e₁e₃ + 3e₂e₃
```
However, conversion can fail if the target type cannot contain the result:
```
julia> convert(KVector{2,VGA(3)}, q)
ERROR: InexactError: convert(KVector{2, VGA(3)}, EvenCliffordNumber{VGA(3), Int64}(0, 1, 2, 3))
...
```

Promotion attempts to preserve the semantics of CliffordNumbers.jl, and therefore prefers to return
an `AbstractCliffordNumber{VGA(3)}`.

```julia-repl
julia> promote_type(EvenCliffordNumber{VGA(3),Int}, QuaternionF64)
EvenCliffordNumber{VGA(3), Float64, 4}

julia> promote_type(KVector{2,VGA(3),Int}, QuaternionF64)
EvenCliffordNumber{VGA(3), Float64, 4}

julia> promote_type(KVector{1,VGA(3),Int}, QuaternionF64)
CliffordNumber{VGA(3), Float64, 8}
```

## [ChainRulesCore.jl]

[ChainRulesCore.jl] defines the `rrule`/`frule` interface that reverse-mode automatic
differentiation engines (Zygote.jl, Diffractor.jl, and others) build on. Loading it alongside
CliffordNumbers.jl registers reverse-mode rules so that scalar losses built from Clifford-number
operations can be differentiated.

Rules are provided for the geometric product, scalar multiplication and division, the grade
automorphisms (`reverse`, `adjoint`, `conj`, `grade_involution`), [`scalar_product`](@ref),
`abs2`, and addition/subtraction. The sandwich product `R * p * R'` (and hence rotor application),
[`versor_inverse`](@ref CliffordNumbers.versor_inverse), and `inv` differentiate by composition
through these rules — no dedicated rule is needed.

The cotangent (gradient) of a Clifford number is itself a Clifford number of the same type, and a
`ChainRulesCore.ProjectTo` is registered to project gradients back onto the grade structure of the
primal. The geometric-product pullback is the exact transpose of the bilinear Cayley product,
computed at compile time from the same machinery as the product, so it is correct in every metric
signature (Euclidean, indefinite, and degenerate). In a positive-definite algebra it reduces to the
familiar `Ω̄ -> (Ω̄ * y', x' * Ω̄)`.

```julia-repl
julia> using CliffordNumbers, Zygote

julia> p = KVector{1,VGA(3)}(1.0, 0.0, 0.0); q = KVector{1,VGA(3)}(0.0, 1.0, 0.0);

julia> loss(R) = abs2(R * p * R' - q);   # squared error of rotating p onto q

julia> R = one(EvenCliffordNumber{VGA(3),Float64});

julia> only(Zygote.gradient(loss, R))    # gradient is an even Clifford number
4-element EvenCliffordNumber{VGA(3), Float64}:
...
```

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

[Quaternions.jl]:   https://github.com/JuliaGeometry/Quaternions.jl
[Unitful.jl]: https://github.com/PainterQubits/Unitful.jl
[ChainRulesCore.jl]: https://github.com/JuliaDiff/ChainRulesCore.jl
