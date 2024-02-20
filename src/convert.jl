import Base.convert

convert(T::Type{<:AbstractCliffordNumber{Q}}, x::AbstractCliffordNumber{Q}) where Q = T(x)::T

function convert(S::Type{<:Real}, x::CliffordNumber)
    isscalar(x) || throw(InexactError(:convert, S, x))
    return S(x[BitIndices(x)[1]])
end
