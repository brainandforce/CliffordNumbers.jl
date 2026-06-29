#---Determinant, trace, and characteristic polynomial----------------------------------------------#
# The determinant, trace, and characteristic polynomial of a multivector `x` are those of its
# left-multiplication operator L_x : y ↦ x*y, a linear map on the 2ⁿ-dimensional algebra. They are
# intrinsic scalar invariants of `x`: in particular `x` is invertible iff `det(x) ≠ 0`, which a
# universal multivector inverse can build on.
#
# Everything is computed matrix-free by the Faddeev–LeVerrier recursion on L_x. Each step is one
# geometric product and one trace, with no 2ⁿ×2ⁿ matrix: the regular representation L is an algebra
# homomorphism, so the recurrence's matrices are themselves L_m for algebra elements `m`, and
# tr(L_y) = 2ⁿ · scalar(y) because only the scalar part of `y` lands on the diagonal. The result
# matches `det`/`tr` of the dense 2ⁿ×2ⁿ matrix in every signature (Euclidean, indefinite, degenerate)
# by construction.

"""
    CliffordNumbers.tr(x::AbstractCliffordNumber{Q,T}) -> T

Returns the trace of the left-multiplication operator `L_x : y ↦ x*y`, equal to `2ⁿ · scalar(x)`
where `n = dimension(Q)`.

When `LinearAlgebra` is loaded, `LinearAlgebra.tr` calls this method.
"""
tr(x::AbstractCliffordNumber) = blade_count(signature(x)) * scalar(x)

# Faddeev–LeVerrier coefficients are integers for an integer carrier; `÷` keeps them in `T` instead
# of widening to a float, and is exact because the algorithm guarantees divisibility.
_leverrier_div(num::Integer, k::Integer) = num ÷ k
_leverrier_div(num, k::Integer) = num / k

"""
    CliffordNumbers.det(x::AbstractCliffordNumber{Q,T}) -> T

Returns the determinant of the left-multiplication operator `L_x : y ↦ x*y`, an intrinsic scalar
invariant of `x`. It satisfies `det(λ) = λ^(2ⁿ)` for a scalar `λ`, `det(x*y) = det(x)*det(y)`, and
`x` is invertible exactly when `det(x) ≠ 0` (so `det(1 + e₁) == 0` in `VGA(2)`).

For the low-dimensional algebras (`2ⁿ ≤ 32`, i.e. `n ≤ 5`) this is the matrix-free Faddeev–LeVerrier
recursion, which is allocation-free on the `Float64` carrier. Above that the generated dense product
becomes too large to compile, so the determinant is instead reduced on the explicit `2ⁿ×2ⁿ` operator
matrix, built from runtime blade arithmetic (no generated kernel), which compiles in O(1) at any
dimension. When `LinearAlgebra` is loaded, `LinearAlgebra.det` calls this method.
"""
det(x::AbstractCliffordNumber{Q}) where Q = _det(Val(blade_count(Q) <= 32), x)

# Low-dimensional algebras: matrix-free Faddeev–LeVerrier on L_x, isbits and allocation-free.
function _det(::Val{true}, x::AbstractCliffordNumber{Q,T}) where {Q,T}
    N = blade_count(Q)
    xd = CliffordNumber{Q,T}(x)
    o = one(xd)
    m = o
    c = -tr(xd)                             # c₁ = -tr(L_x)
    for k in 2:N
        m = muladd(c, o, xd * m)           # Mₖ = L_x Mₖ₋₁ + cₖ₋₁ I
        c = _leverrier_div(-tr(xd * m), k) # cₖ = -tr(L_x Mₖ) / k
    end
    return ifelse(isodd(N), -c, c)         # det = (-1)ᴺ c_N
end

# High-dimensional algebras: reduce the explicit left-multiplication matrix. Avoids compiling the
# 2ⁿ×2ⁿ-term generated geometric product, which is prohibitive past n = 5.
_det(::Val{false}, x::AbstractCliffordNumber) = _bareiss_det(_leftmul_matrix(x))

# Matrix of L_x : y ↦ x*y in the blade basis, built without the generated product. The coefficient of
# eᵢ in x*eⱼ comes solely from the eₖ component of x with k = i ⊻ j, so L[i,j] = x_{i⊻j} · sign(eₖeⱼ).
function _leftmul_matrix(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    N = blade_count(Q)
    M = Matrix{T}(undef, N, N)
    for j in 0:(N - 1)
        bj = BladeIndex{Q}(UInt(j))
        for i in 0:(N - 1)
            bk = BladeIndex{Q}(UInt(xor(i, j)))
            M[i + 1, j + 1] = x[bk] * sign_of_mult(bk, bj)
        end
    end
    return M
end

# Bareiss fraction-free elimination: exact for integer carriers (every division is exact) and for
# the rational/float fields, keeping the result in the scalar type.
_exact_div(a::Integer, b::Integer) = a ÷ b
_exact_div(a, b) = a / b

function _bareiss_det(M::Matrix{T}) where T
    n = size(M, 1)
    s = 1
    prev = one(T)
    for k in 1:(n - 1)
        if iszero(M[k, k])
            r = findnext(!iszero, view(M, :, k), k + 1)
            r === nothing && return zero(T)
            for c in 1:n
                M[k, c], M[r, c] = M[r, c], M[k, c]
            end
            s = -s
        end
        for j in (k + 1):n, i in (k + 1):n
            M[i, j] = _exact_div(M[i, j] * M[k, k] - M[i, k] * M[k, j], prev)
        end
        prev = M[k, k]
    end
    return s * M[n, n]
end

"""
    CliffordNumbers.charpoly_coeffs(x::AbstractCliffordNumber{Q,T}) -> NTuple{2ⁿ,T}

Returns the coefficients `(c₁, …, c_N)` of the characteristic polynomial of `L_x`,
`p(λ) = λᴺ + c₁ λᴺ⁻¹ + ⋯ + c_N` with `N = 2ⁿ`. The trace is `-c₁` and the determinant is
`(-1)ᴺ c_N`. By the Cayley–Hamilton theorem `x` satisfies its own characteristic polynomial:
`xᴺ + c₁ xᴺ⁻¹ + ⋯ + c_N == 0`.
"""
@generated function charpoly_coeffs(x::AbstractCliffordNumber{Q,T}) where {Q,T}
    N = blade_count(Q)
    body = quote
        xd = CliffordNumber{Q,T}(x)
        o = one(xd)
        m = o
        c_1 = -tr(xd)
    end
    names = Any[:c_1]
    for k in 2:N
        ck = Symbol(:c_, k)
        ckm1 = Symbol(:c_, k - 1)
        push!(body.args, quote
            m = muladd($ckm1, o, xd * m)
            $ck = _leverrier_div(-tr(xd * m), $k)
        end)
        push!(names, ck)
    end
    push!(body.args, :(return tuple($(names...))))
    return body
end
