import Base.promote_rule

# Generic promote rule for Clifford numbers with different element types
function promote_rule(
    ::Type{<:CliffordNumber{Q,T1}},
    ::Type{<:CliffordNumber{Q,T2}}
) where {Q,T1,T2}
    T = promote_type(T1,T2)
    return CliffordNumber{Q,T,elements(Q)}
end

# Promote rule when using real numbers
function promote_rule(
    ::Type{<:CliffordNumber{Q,T1}},
    T2::Type{<:Real}
) where {Q,T1}
    T = promote_type(T1,T2)
    return CliffordNumber{Q,T,elements(Q)}
end

# Promote rule for complex numbers with real Clifford numbers
function promote_rule(
    ::Type{<:CliffordNumber{Q,T1}},
    ::Type{<:Complex{T2}}
) where {Q,T1<:Real,T2}
    T = promote_type(T1,T2)
    return CliffordNumber{Q,T,elements(Q)}
end

# For abstract types where no element type of a Clifford number is specified
promote_rule(::Type{CliffordNumber{Q}}, S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S
promote_rule(::Type{CliffordNumber{Q}}, ::Type{T}) where {Q,T<:Real} = CliffordNumber{Q,T}

# TODO: may need macros to handle the case of a CliffordNumber{Cl} being promoted with a Complex
# In that case, promotion depends on whether the pseudoscalar behaves like the imaginary unit
