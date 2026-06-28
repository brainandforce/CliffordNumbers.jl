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

"""
    CliffordNumbers._complement_expr(::Type{T}, blade_complement) where T<:AbstractCliffordNumber

Builds the body of a compile-time complement of `x::T`. `blade_complement` is the `BladeIndex`-level
complement applied to each output blade index, mirroring the definitions

    left_complement(x)  = C(x[right_complement.(BladeIndices(C))])
    right_complement(x) = C(x[left_complement.(BladeIndices(C))])

where `C = complement_type(T)`. Resolving each output coefficient to a tuple position and sign at
compile time avoids the runtime `to_index` blade search the generic indexing path incurs for
`KVector` (and the dynamic dispatch that path emits, which does not lower on GPU).
"""
function _complement_expr(::Type{T}, blade_complement) where T<:AbstractCliffordNumber
    C = complement_type(T)
    grades = nonzero_grades(T)
    terms = Any[]
    # Iterate output blades in C's storage order; `b` is the source blade in `x`
    for c in BladeIndices(C)
        b = blade_complement(c)
        if grade(b) in grades
            # x[b] == flipsign(x.data[to_index(T, b)], b); resolve both at compile time
            push!(terms, :(flipsign(Tuple(x)[$(to_index(T, b))], $(sign(b)))))
        else
            push!(terms, :(zero(scalar_type(x))))
        end
    end
    return :($C(tuple($(terms...))))
end

@generated left_complement(x::AbstractCliffordNumber) = _complement_expr(x, right_complement)
@generated right_complement(x::AbstractCliffordNumber) = _complement_expr(x, left_complement)
