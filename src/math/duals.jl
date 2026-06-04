#---Sign changing operations-----------------------------------------------------------------------#

"""
    CliffordNumbers._sign_automorphism(x::AbstractCliffordNumber, ::Val{F})

Applies the sign-changing grade automorphism named by the `Symbol` `F` (`:adjoint`, `:reverse`,
`:conj`, or `:grade_involution`) to `x`.

These automorphisms act blade-by-blade: each coefficient maps to a fixed storage position with a
fixed sign. The generic indexing form `T(x[F.(BladeIndices(T))])` resolves those positions at
runtime through `to_index`, which does not lower on the GPU (it allocates and dispatches
dynamically). This generated implementation resolves every position and sign at compile time and
emits pure tuple arithmetic, which is allocation-free and GPU-safe.
"""
@generated function _sign_automorphism(x::T, ::Val{F}) where {T<:AbstractCliffordNumber,F}
    f = getfield(@__MODULE__, F)
    inds = BladeIndices(T)
    elements = map(1:length(inds)) do i
        b = f(inds[i])
        return :($(sign(b)) * Tuple(x)[$(to_index(T, b))])
    end
    return :(T(($(elements...),)))
end

for f in (:adjoint, :conj, :grade_involution, :reverse)
    @eval $f(x::AbstractCliffordNumber) = _sign_automorphism(x, Val($(QuoteNode(f))))
end

# Faster implementations for KVector that don't require indexing. For real scalars `reverse` and
# `adjoint` apply the same per-grade sign, so they share an implementation.
adjoint(x::KVector{K}) where K = ifelse(iszero(K & 2), x, -x)
reverse(x::KVector{K}) where K = ifelse(iszero(K & 2), x, -x)
grade_involution(x::KVector{K}) where K = ifelse(iseven(K), x, -x)
conj(x::KVector{K}) where K = ifelse(iszero((K + 1) & 2), x, -x)

function (::CyclicGradeNegation{N,I,R})(x::AbstractCliffordNumber) where {N,I,R}
    x = ifelse(N, x, -x)
    x = ifelse(I, x, grade_involution(x))
    x = ifelse(R, x, reverse(x))
    return x
end

function left_complement(x::T) where T<:AbstractCliffordNumber
    C = complement_type(T)
    return C(x[right_complement.(BladeIndices(C))])
end

function right_complement(x::T) where T<:AbstractCliffordNumber
    C = complement_type(T)
    return C(x[left_complement.(BladeIndices(C))])
end
