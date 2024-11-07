module CliffordNumbersQuaternionsExt

using CliffordNumbers
using Quaternions

CliffordNumbers.scalar_type(::Type{Quaternion{T}}) where T = T

#---Conversion methods-----------------------------------------------------------------------------#

# To AbstractCliffordNumber
(::Type{C})(q::Quaternion) where C<:EvenCliffordNumber{VGA(3)} = C(q.s, q.v1, q.v2, q.v3)
AbstractCliffordNumber(q::Quaternion) = EvenCliffordNumber{VGA(3)}(q)

function (::Type{C})(q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return C(EvenCliffordNumber{VGA(3)}(q))
end

function Base.convert(::Type{C}, q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return convert(C, EvenCliffordNumber{VGA(3)}(q))
end

# To Quaternion
"""
    Quaternion(c::AbstractCliffordNumber{VGA(3)})
    Quaternion{T}(c::AbstractCliffordNumber{VGA(3)})

Constructs a quaternion from an element of the algebra of physical space, the 3D geometric algebra
with a positive-definite signature whose even subalgebra is isomorphic to the quaternion algebra ℍ.
"""
(::Type{H})(c::EvenCliffordNumber{VGA(3)}) where H<:Quaternion = H(Tuple(c)...)

function (::Type{H})(c::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
    return H(Tuple(EvenCliffordNumber{VGA(3)}(c))...)
end

function Base.convert(::Type{H}, c::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
    return H(Tuple(convert(EvenCliffordNumber{VGA(3)}, c))...)
end

# Type promotion
function Base.promote_rule(
    ::Type{C}, 
    ::Type{Quaternion{T}}
) where {S,C<:AbstractCliffordNumber{VGA(3),S},T}
    return promote_type(C, EvenCliffordNumber{VGA(3),T})
end

#---Arithmetic operations--------------------------------------------------------------------------#

end
