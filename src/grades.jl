#---Tools for determining which grades of a Clifford number are nonzero----------------------------#
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

"""
    has_grades_of(S::Type{<:AbstractCliffordNumber}, T::Type{<:AbstractCliffordNumber}) -> Bool
    has_grades_of(x::AbstractCliffordNumber, y::AbstractCliffordNumber) -> Bool

Returns `true` if the grades represented in S are also represented in T; `false` otherwise.
"""
function has_grades_of(
    ::Type{S},
    ::Type{T}
) where {S<:AbstractCliffordNumber,T<:AbstractCliffordNumber}
    return all(x in nonzero_grades(T) for x in nonzero_grades(S))
end

function has_grades_of(x::AbstractCliffordNumber, T::Type{<:AbstractCliffordNumber})
    return has_grades_of(typeof(x), T)
end

function has_grades_of(T::Type{<:AbstractCliffordNumber}, x::AbstractCliffordNumber)
    return has_grades_of(T, typeof(x))
end

function has_grades_of(x::AbstractCliffordNumber, y::AbstractCliffordNumber)
    return has_grades_of(typeof(x), typeof(y))
end
