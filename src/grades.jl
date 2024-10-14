#---Tools for determining which grades of a Clifford number are nonzero----------------------------#
"""
    nonzero_grades(::Type{<:AbstractCliffordNumber{Q}}) -> AbstractVector{Int}
    nonzero_grades(::Number)

Returns an `AbstractVector{Int}` whose elements are all nonzero grades of a Clifford number or type.
Any subtype `T` of `AbstractCliffordNumber` should define this method for `Type{T}`; it is
automatically implemented for instance arguments.

For arguments that are `Real` or `Complex`, this function returns `0:0`.

# Examples

```julia-repl
julia> CliffordNumbers.nonzero_grades(CliffordNumber{VGA(3)})
0:3

julia> CliffordNumbers.nonzero_grades(EvenCliffordNumber{VGA(3)})
0:2:2

julia> CliffordNumbers.nonzero_grades(KVector{2,VGA(3)})
2:2

julia> nonzero_grades(Float64)
0:0
```
"""
nonzero_grades(x::Number) = nonzero_grades(typeof(x))
nonzero_grades(::Type{<:BaseNumber}) = 0:0

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
