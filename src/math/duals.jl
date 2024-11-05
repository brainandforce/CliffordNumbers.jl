#---Sign changing operations-----------------------------------------------------------------------#

for f in (:adjoint, :grade_involution, :conj)
    @eval $f(x::T) where T<:AbstractCliffordNumber = x[$f.(BitIndices(T))]
end

reverse(x::AbstractCliffordNumber) = adjoint(x)

# Faster implementations for KVector that don't require indexing
adjoint(x::KVector{K}) where K = ifelse(iszero(K & 2), x, -x)
grade_involution(x::KVector{K}) where K = ifelse(iseven(K), x, -x)
conj(x::KVector{K}) where K = ifelse(iszero((K + 1) & 2), x, -x)

left_complement(x::AbstractCliffordNumber) = x[right_complement.(BitIndices(complement_type(x)))]
right_complement(x::AbstractCliffordNumber) = x[left_complement.(BitIndices(complement_type(x)))]
