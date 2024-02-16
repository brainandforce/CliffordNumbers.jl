import Base.convert

convert(S::Type{<:CliffordNumber{Q,T}}, x::CliffordNumber{Q}) where {Q,T} = S(x.data)::S

function convert(S::Type{<:Real}, x::CliffordNumber)
    isscalar(x) || throw(InexactError(:convert, S, x))
    return S(x[BitIndices(x)[1]])
end
