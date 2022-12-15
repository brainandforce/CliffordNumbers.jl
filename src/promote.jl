import Base.promote_rule

# Generic promote rule for Clifford numbers with different element types
function promote_rule(
    ::Type{<:CliffordNumber{Cl,T1}},
    ::Type{<:CliffordNumber{Cl,T2}}
) where {Cl,T1,T2}
    T = promote_type(T1,T2)
    return CliffordNumber{Cl,T,elements(Cl)}
end

# Promote rule when using real numbers
function promote_rule(
    ::Type{<:CliffordNumber{Cl,T1}},
    T2::Type{<:Real}
) where {Cl,T1}
    T = promote_type(T1,T2)
    return CliffordNumber{Cl,T,elements(Cl)}
end

# Promote rule for complex numbers with real Clifford numbers
function promote_rule(
    ::Type{<:CliffordNumber{Cl,T1}},
    ::Type{<:Complex{T2}}
) where {Cl,T1<:Real,T2}
    T = promote_type(T1,T2)
    return CliffordNumber{Cl,T,elements(Cl)}
end

# For abstract types where no element type of a Clifford number is specified
promote_rule(::Type{CliffordNumber{Cl}}, S::Type{<:CliffordNumber{Cl,T}}) where {Cl,T} = S
promote_rule(::Type{CliffordNumber{Cl}}, ::Type{T}) where {Cl,T<:Real} = CliffordNumber{Cl,T}

# TODO: may need macros to handle the case of a CliffordNumber{Cl} being promoted with a Complex
# In that case, promotion depends on whether the pseudoscalar behaves like the imaginary unit
