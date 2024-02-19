#---Tools for determining which grades of a Clifford number are nonzero----------------------------#
"""
    RepresentedGrades{C<:AbstractCliffordNumber} <: AbstractVector{Bool}

A vector describing the grades explicitly represented by some Clifford number type `C`.

Indexing of a `RepresentedGrades{Q}` is zero-based, with indices corresponding to all grades from
0 (scalar) to `dimension(Q)` (pseudoscalar), with `true` values indicating nonzero grades, and
`false` indicating that the grade is zero.

# Implementation

The function `nonzero_grades(::Type{T})` should be implemented for types `T` descending from 
`AbstractCliffordNumber`.
"""
struct RepresentedGrades{C<:AbstractCliffordNumber} <: AbstractVector{Bool}
end

RepresentedGrades(x::AbstractCliffordNumber) = RepresentedGrades{typeof(x)}()
RepresentedGrades(C::Type{<:AbstractCliffordNumber}) = RepresentedGrades(zero(C))

Base.size(::RepresentedGrades{C}) where C = tuple(dimension(QuadraticForm(C)) + 1)
Base.axes(::RepresentedGrades{C}) where C = tuple(0:dimension(QuadraticForm(C)))

"""
    nonzero_grades(::Type{<:AbstractCliffordNumber})
    nonzero_grades(::AbstractCliffordNumber)

A function returning an indexable object representing all nonzero grades of a Clifford number
representation.

This function is used to define the indexing of `RepresentedGrades`, and should be defined for any
subtypes of `AbstractCliffordNumber`.

# Examples

```julia-repl
julia> CliffordNumbers.nonzero_grades(CliffordNumber{APS})
0:3

julia> CliffordNumbers.nonzero_grades(KVector{2,APS})
2:2
```
"""
nonzero_grades(x::Number) = nonzero_grades(typeof(x))
# TODO: define for Complex?
nonzero_grades(::Type{<:Real}) = 0:0

function Base.getindex(r::RepresentedGrades{C}, i::Integer) where C
    @boundscheck checkbounds(r, i)
    return i in nonzero_grades(C)
end
