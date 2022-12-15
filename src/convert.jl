import Base.convert

convert(S::Type{<:CliffordNumber{Cl,T}}, m::CliffordNumber{Cl}) where {Cl,T} = S(m.data)::S
