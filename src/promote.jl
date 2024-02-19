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
promote_rule(C::Type{<:AbstractCliffordNumber}, N::Type{<:BaseNumber}) = similar_type(C, N)
