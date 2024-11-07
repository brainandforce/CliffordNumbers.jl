#---Sign changing operations-----------------------------------------------------------------------#

for f in (:adjoint, :conj, :grade_involution, :reverse)
    @eval $f(x::T) where T<:AbstractCliffordNumber = T(x[$f.(BladeIndices(T))])
end

# Faster implementations for KVector that don't require indexing
adjoint(x::KVector{K}) where K = ifelse(iszero(K & 2), x, -x)
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
