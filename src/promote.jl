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
function promote_rule(
    S::Type{<:AbstractCliffordNumber{Q}},
    T::Type{<:AbstractCliffordNumber{Q}}
) where Q
    return CliffordNumber{Q,promote_numeric_type(S,T),elements(Q)}
end

# Promote rule for BaseNumber types (real, complex)
# Note that complex numbers aren't automatically treated as pseudoscalars
# (this only works in some dimensions...)
function promote_rule(::Type{<:CliffordNumber{Q,T}}, ::Type{N}) where {Q,T,N<:BaseNumber}
    return CliffordNumber{Q,promote_type(T,N),elements(Q)}
end

function promote_rule(::Type{CliffordNumber{Q}}, ::Type{N}) where {Q,N<:BaseNumber}
    return CliffordNumber{Q,N,elements(Q)}
end

#---Promotion rules for various representations----------------------------------------------------#

function promote_rule(S::Type{<:KVector{K,Q}}, T::Type{<:KVector{K,Q}}) where {K,Q}
    return KVector{K,Q,promote_numeric_type(S,T),binomial(dimension(Q),K)}
end

function promote_rule(S::Type{<:KVector{K1,Q}}, T::Type{<:KVector{K2,Q}}) where {K1,K2,Q}
    # TODO: logic for even/odd multivector types
    return CliffordNumber{Q,promote_numeric_type(S,T),elements(Q)}
end

function promote_rule(::Type{<:KVector{<:Any,Q,T}}, ::Type{N}) where {Q,T,N<:BaseNumber}
    return CliffordNumber{Q,promote_type(T,N),elements(Q)}
end

# 0-vectors are scalars, but keep CliffordNumbers semantics.
function promote_rule(::Type{<:KVector{0,Q,T}}, ::Type{N}) where {Q,T,N<:BaseNumber}
    return KVector{0,Q,promote_type(T,N),1}
end
