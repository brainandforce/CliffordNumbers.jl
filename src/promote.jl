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
        return :(Z2CliffordNumber{isodd(K1),Q,promote_type(T1,T2),div(elements(Q),2)})
    end
end

function promote_rule(C::Type{<:KVector{K,Q,T}}, ::Type{N}) where{K,Q,T,N<:BaseNumber}
    return promote_type(C, KVector{0,Q,N,1})
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

function promote_rule(C::Type{<:Z2CliffordNumber{P,Q,T}}, ::Type{N}) where {P,Q,T,N<:BaseNumber}
    return promote_type(C, KVector{0,Q,N,1})
end

#= TODO:
    This needs to be defined because of how Base.promote_type is implemented and how the Default
    type promotion for CliffordNumber subtypes is defined.
=#
promote_rule(C1::Type{<:KVector}, C2::Type{<:Z2CliffordNumber}) = promote_rule(C2, C1)

#---Promoting just the scalar types of an AbstractCliffordNumber-----------------------------------#
"""
    scalar_promote(x::AbstractCliffordNumber, y::AbstractCliffordNumber)

Promotes the scalar types of `x` and `y` to a common type. This does not increase the number of
represented grades of either `x` or `y`.
"""
scalar_promote() = ()
scalar_promote(x::Number) = tuple(x)

function scalar_promote(x::Number, y::Number)
    @inline
    T = promote_numeric_type(x, y)
    return (scalar_convert(T, x), scalar_convert(T, y))
end

scalar_promote(x::T, y::T) where T = (x,y)
scalar_promote(x::T, y::AbstractCliffordNumber{<:Any,T}) where T = (x,y)
scalar_promote(x::AbstractCliffordNumber{<:Any,T}, y::T) where T = (x,y)

function scalar_promote(
    x::AbstractCliffordNumber{<:Any,T},
    y::AbstractCliffordNumber{<:Any,T}
) where T
    return (x,y)
end

function scalar_promote(x::Number, y::Number, zs::Number...)
    @inline
    T = promote_numeric_type(x, y, zs...)
    return scalar_convert.(T, (x, y, zs...))
end

#---Widening types---------------------------------------------------------------------------------#
"""
    widen(C::Type{<:AbstractCliffordNumber})
    widen(x::AbstractCliffordNumber)

Construct a new type whose scalar type is widened. This behavior matches that of
`widen(C::Type{Complex{T}})`, which results in widening of its scalar type `T`.

For obtaining a representation of a Clifford number with an increased number of nonzero grades,
use `widen_grade(T)`.
"""
Base.widen(C::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = similar_type(C, widen(T))

"""
    widen_grade(C::Type{<:AbstractCliffordNumber})
    widen_grade(x::AbstractCliffordNumber)

For type arguments, construct the next largest type that can hold all of the grades of `C`.
`KVector{K,Q,T}` widens to `EvenCliffordNumber{Q,T}` or `OddCliffordNumber{Q,T}`, and
`EvenCliffordNumber{Q,T}` and `OddCliffordNumber{Q,T}` widen to `CliffordNumber{Q,T}`, which is the
widest type.

For `AbstractCliffordNumber` arguments, the argument is converted to the result of
`widen_grade(typeof(x))`.

For widening the scalar type of an `AbstractCliffordNumber`, use `Base.widen(T)`.
"""
widen_grade(x::AbstractCliffordNumber) = convert(widen_grade(typeof(x)), x)

widen_grade(C::Type{<:CliffordNumber}) = C

widen_grade(::Type{Z2CliffordNumber}) = CliffordNumber
widen_grade(::Type{Z2CliffordNumber{P}}) where {P} = CliffordNumber
widen_grade(::Type{Z2CliffordNumber{P,Q}}) where {P,Q} = CliffordNumber{Q}
widen_grade(::Type{Z2CliffordNumber{P,Q,T}}) where {P,Q,T} = CliffordNumber{Q,T}
widen_grade(::Type{Z2CliffordNumber{P,Q,T,L}}) where {P,Q,T,L} = CliffordNumber{Q,T,elements(Q)}

widen_grade(::Type{KVector}) = CliffordNumber
widen_grade(::Type{KVector{K}}) where {K} = Z2CliffordNumber{isodd(K)}
widen_grade(::Type{KVector{K,Q}}) where {K,Q} = Z2CliffordNumber{isodd(K),Q}
widen_grade(::Type{KVector{K,Q,T}}) where {K,Q,T} = Z2CliffordNumber{isodd(K),Q,T}

function widen_grade(::Type{KVector{K,Q,T,L}}) where {K,Q,T,L}
    return Z2CliffordNumber{isodd(K),Q,T,div(elements(Q), 2)}
end
