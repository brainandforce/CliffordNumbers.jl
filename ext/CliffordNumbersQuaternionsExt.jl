module CliffordNumbersQuaternionsExt

using CliffordNumbers
using Quaternions

import Base: convert, promote_rule
import Base: *

import CliffordNumbers: AbstractCliffordNumber, EvenCliffordNumber
import CliffordNumbers: ∧, ∨, ⨼, ⨽, ×, ⨰, dot, scalar_product

CliffordNumbers.scalar_type(::Type{Quaternion{T}}) where T = T

#---Conversion methods-----------------------------------------------------------------------------#

# To AbstractCliffordNumber
(::Type{C})(q::Quaternion) where C<:EvenCliffordNumber{VGA(3)} = C(q.s, q.v1, q.v2, q.v3)
AbstractCliffordNumber(q::Quaternion) = EvenCliffordNumber{VGA(3)}(q)

function (::Type{C})(q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return C(EvenCliffordNumber{VGA(3)}(q))
end

function convert(::Type{C}, q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return convert(C, EvenCliffordNumber{VGA(3)}(q))
end

# To Quaternion
"""
    Quaternion(c::AbstractCliffordNumber{VGA(3)})
    Quaternion{T}(c::AbstractCliffordNumber{VGA(3)})

Constructs a quaternion from an element of the algebra of physical space, the 3D geometric algebra
with a positive-definite signature whose even subalgebra is isomorphic to the quaternion algebra ℍ.
Any odd-grade coefficients of `c` are lost.

If loss of odd-grade coefficients should throw an error, use `convert(Quaternion, c)` or
`convert(Quaternion{T}, c)` instead of the constructor.
"""
(::Type{H})(c::EvenCliffordNumber{VGA(3)}) where H<:Quaternion = H(Tuple(c)...)

function (::Type{H})(c::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
    return H(EvenCliffordNumber{VGA(3)}(c))
end

function convert(::Type{H}, c::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
    # This fails if Quaternion would drop represented grades of c
    return H(convert(EvenCliffordNumber{VGA(3)}, c))
end

# Type promotion
function promote_rule(
    ::Type{C}, 
    ::Type{Quaternion{T}}
) where {S,C<:AbstractCliffordNumber{VGA(3),S},T}
    return promote_type(C, EvenCliffordNumber{VGA(3),T})
end

#---Arithmetic operations--------------------------------------------------------------------------#

for op in (:*, :∧, :∨, :⨼, :⨽, :×, :⨰, :dot, :scalar_product)
    @eval begin
        function $op(q::Quaternion, x::AbstractCliffordNumber{Q}) where Q
            return $op(convert(AbstractCliffordNumber{Q}, q), x)
        end
        function $op(x::AbstractCliffordNumber{Q}, q::Quaternion) where Q
            return $op(x, convert(AbstractCliffordNumber{Q}, q))
        end
    end
end

end
