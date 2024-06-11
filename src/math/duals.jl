#---Sign changing operations-----------------------------------------------------------------------#

for f in (:adjoint, :grade_involution, :conj)
    @eval $f(x::T) where T<:AbstractCliffordNumber = x[$f.(BitIndices(T))]
end

reverse(x::AbstractCliffordNumber) = adjoint(x)

# Faster implementations for KVector that don't require indexing
adjoint(x::T) where T<:KVector = T(x.data .* Int8(-1)^!iszero(grade(x) & 2))
grade_involution(x::T) where T<:KVector = T(x.data .* Int8(-1)^isodd(grade(x)))
conj(x::T) where T<:KVector = T(x.data .* Int8(-1)^!iszero((grade(x) + 1) & 2))

left_complement(x::AbstractCliffordNumber) = x[right_complement.(BitIndices(complement_type(x)))]
right_complement(x::AbstractCliffordNumber) = x[left_complement.(BitIndices(complement_type(x)))]
