import Base.convert

# Default conversion should check for exact representability
function convert(::Type{T}, x::AbstractCliffordNumber{Q}) where {Q,T<:AbstractCliffordNumber{Q}}
    has_grades_of(x, T) && return T(x)::T
    throw(InexactError(:convert, T, x))
end

function convert(::Type{T}, x::AbstractCliffordNumber) where T<:BaseNumber
    return isscalar(x) ? T(x[BitIndex{QuadraticForm(x)}()]) : throw(InexactError(:convert, T, x))
end
