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
function promote_rule(S::Type{<:CliffordNumber{Q}}, T::Type{<:CliffordNumber{Q}}) where Q
    return CliffordNumber{Q,promote_eltypeof(S,T),elements(Q)}
end

# Promote rule when using real numbers
function promote_rule(S::Type{<:CliffordNumber{Q,<:Real}}, T::Type{<:Real}) where Q
    return CliffordNumber{Q,promote_eltypeof(S,T),elements(Q)}
end

# Promote rule for complex numbers with real Clifford numbers
function promote_rule(S::Type{<:CliffordNumber{Q,<:Real}}, T::Type{<:Complex{R}}) where {Q,R}
    return CliffordNumber{Q,promote_eltypeof(S,R),elements(Q)}
end

# For abstract types where no element type of a Clifford number is specified
promote_rule(::Type{CliffordNumber{Q}}, S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S
promote_rule(::Type{CliffordNumber{Q}}, ::Type{T}) where {Q,T<:Real} = CliffordNumber{Q,T}

# TODO: may need macros to handle the case of a CliffordNumber{Cl} being promoted with a Complex
# In that case, promotion depends on whether the pseudoscalar behaves like the imaginary unit
