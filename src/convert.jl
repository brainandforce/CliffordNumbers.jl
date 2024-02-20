import Base.convert

# Default conversion should check for exact representability
function convert(T::Type{<:AbstractCliffordNumber{Q}}, x::AbstractCliffordNumber{Q}) where Q
    has_grades_of(x, T) && return T(x)
    throw(InexactError(:convert, T, x))
end

function convert(S::Type{<:Real}, x::CliffordNumber)
    isscalar(x) || throw(InexactError(:convert, S, x))
    return S(x[BitIndices(x)[1]])
end
