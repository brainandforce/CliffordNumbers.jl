import Base.convert

convert(S::Type{<:CliffordNumber{Q,T}}, m::CliffordNumber{Q}) where {Q,T} = S(m.data)::S
