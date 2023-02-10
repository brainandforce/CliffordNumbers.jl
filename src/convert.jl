import Base.convert

convert(S::Type{<:CliffordNumber{Q,T}}, m::CliffordNumber{Q}) where {Q,T} = S(m.data)::S

function convert(S::Type{<:Real}, m::CliffordNumber)
    isscalar(m) || throw(InexactError(:convert, S, m))
    return S(m[0])
end
