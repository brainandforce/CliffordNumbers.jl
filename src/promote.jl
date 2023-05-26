import Base.promote_rule

# This set of functions is defined in Julia Base (abstractarray.jl), but rather than importing it
# we're duplicating it (this part of the Julia API may be unstable...)
promote_eltype() = Base.Bottom
promote_eltype(x, y...) = promote_type(eltype(x), promote_eltype(y...))

eltypeof(x) = typeof(x)
eltypeof(x::AbstractCliffordNumber) = eltype(x)

promote_eltypeof() = Base.Bottom
promote_eltypeof(x, y...) = promote_type(eltypeof(x), promote_eltypeof(y...))

# Generic promote rule for Clifford numbers with different element types
function promote_rule(
    S::Type{<:AbstractCliffordNumber{Q}},
    T::Type{<:AbstractCliffordNumber{Q}}
) where Q
    return CliffordNumber{Q,promote_eltype(S,T),elements(Q)}
end

# Promote rule for BaseNumber types (real, complex)
# Note that complex numbers aren't automatically treated as pseudoscalars
# (this only works in some dimensions...)
function promote_rule(S::Type{<:AbstractCliffordNumber{Q}}, T::Type{<:BaseNumber}) where Q
    return CliffordNumber{Q,promote_type(eltype(S),T),elements(Q)}
end
