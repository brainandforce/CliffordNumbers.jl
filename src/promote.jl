import Base.promote_rule

# Inspired by Base.promote_eltype in Julia Base
"""
    promote_numeric_type(x, y)

Calls `promote_type()` for the result of `numeric_type()` called on all arguments. For incompletely
specified types, the result of `numeric_type()` is replaced with `Bool`, which always promotes to
any larger numeric type.
"""
promote_numeric_type() = Union{}

function promote_numeric_type(x, y...)
    T = numeric_type(x)
    return promote_type(ifelse(T >: BaseNumber, Bool, T), promote_numeric_type(y...))
end

# Generic promote rule for Clifford numbers with different element types
#= TODO:
    This definition can screw things up sometimes with promote_type, because it evaluates
    promote_rule in both orders and then promotes the result of both. In the case of KVector and
    Z2CliffordNumber, both orders need to be defined or the result is always CliffordNumber (it 
    should be Z2CliffordNumber if the KVector matches the parity of the Z2CliffordNumber)
=#
function promote_rule(
    ::Type{<:AbstractCliffordNumber{Q,S}},
    ::Type{<:AbstractCliffordNumber{Q,T}}
) where {Q,S,T}
    return CliffordNumber{Q,promote_type(S,T),elements(Q)}
end

# Promote rule for BaseNumber types (real, complex)
# Note that complex numbers aren't automatically treated as pseudoscalars
# (this only works in some dimensions...)
function promote_rule(C::Type{<:AbstractCliffordNumber{Q}}, ::Type{N}) where {Q,N<:BaseNumber}
    T = numeric_type(C)
    return ifelse(T >: BaseNumber, C, similar_type(C, promote_type(T,N)))
end

#---Promotion rules for various representations----------------------------------------------------#

function promote_rule(::Type{<:KVector{K,Q,S}}, ::Type{<:KVector{K,Q,T}}) where {K,Q,S,T}
    return KVector{K,Q,promote_type(S,T),binomial(dimension(Q),K)}
end

@generated function promote_rule(
    ::Type{<:KVector{K1,Q,T1}},
    ::Type{<:KVector{K2,Q,T2}}
) where {K1,K2,Q,T1,T2}
    if xor(iseven(K1), iseven(K2))
        return :(CliffordNumber{Q,promote_type(T1,T2),elements(Q)})
    else
        return :(Z2CliffordNumber{iseven(K1),Q,promote_type(T1,T2),div(elements(Q),2)})
    end
end

function promote_rule(::Type{<:KVector{<:Any,Q,T}}, ::Type{N}) where {Q,T,N<:BaseNumber}
    return CliffordNumber{Q,promote_type(T,N),elements(Q)}
end

# 0-vectors are scalars, but keep CliffordNumbers semantics.
function promote_rule(::Type{<:KVector{0,Q,T}}, ::Type{N}) where {Q,T,N<:BaseNumber}
    return KVector{0,Q,promote_type(T,N),1}
end

function promote_rule(
    ::Type{<:Z2CliffordNumber{P,Q,S}},
    ::Type{<:Z2CliffordNumber{P,Q,T}}
) where {P,Q,S,T}
    return Z2CliffordNumber{P,Q,promote_type(S,T)}
end

@generated function promote_rule(
    ::Type{<:Z2CliffordNumber{P,Q,S}},
    ::Type{<:KVector{K,Q,T}}
) where {P,K,Q,S,T}
    if xor(P, isodd(K))
        return :(CliffordNumber{Q,promote_type(S,T),elements(Q)})
    else
        return :(Z2CliffordNumber{P,Q,promote_type(S,T),div(elements(Q), 2)})
    end
end

#= TODO:
    This needs to be defined because of how Base.promote_type is implemented and how the Default
    type promotion for CliffordNumber subtypes is defined.
=#
promote_rule(C1::Type{<:KVector}, C2::Type{<:Z2CliffordNumber}) = promote_rule(C2, C1)
