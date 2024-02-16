#---Tools for determining which grades of a Clifford number are nonzero----------------------------#
"""
    RepresentedGrades{Q} <: AbstractVector{Bool}

Determines which grades of a Clifford number with quadratic form `Q` are represented by a specific
data type.

The backing `UInt` is constructed by summing the base 2 exponentials of each represented grade. In
this sense, `RepresentedGrades{Q}` is a bit vector with length equal to `dimension(Q) + 1`. All bits
above the maximum grade associated with `Q` are zero.

Indexing of a `RepresentedGrades{Q}` is zero-based, with indices corresponding to all grades from
0 (scalar) to `dimension(Q)` (pseudoscalar).
"""
struct RepresentedGrades{Q} <: AbstractVector{Bool}
    bits::UInt
    RepresentedGrades{Q}(bits::Integer) where Q = new(bits & (2^(dimension(Q) + 1) - 1))
end

Base.size(::RepresentedGrades{Q}) where Q = tuple(dimension(Q) + 1)
Base.axes(::RepresentedGrades{Q}) where Q = tuple(0:dimension(Q))

function Base.getindex(r::RepresentedGrades, i::Integer)
    @boundscheck checkbounds(r, i)
    return !iszero(r.bits & UInt(2)^i)
end

#---Getting represented grades of AbstractCliffordNumber subtypes----------------------------------#

RepresentedGrades(x::AbstractCliffordNumber) = RepresentedGrades(typeof(x))
