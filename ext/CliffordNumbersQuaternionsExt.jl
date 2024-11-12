module CliffordNumbersQuaternionsExt

using CliffordNumbers
using Quaternions

import Base: convert, promote_rule
import Base: *

import CliffordNumbers: AbstractCliffordNumber, EvenCliffordNumber
import CliffordNumbers: ∧, ∨, ⨼, ⨽, ×, ⨰, dot, scalar_product

import Quaternions: slerp

CliffordNumbers.scalar_type(::Type{Quaternion{T}}) where T = T

#---Conversion to Clifford numbers-----------------------------------------------------------------#

# By default, treat a Quaternion{T} like it is an EvenCliffordNumber{VGA(3),T}
(::Type{C})(q::Quaternion) where C<:EvenCliffordNumber{VGA(3)} = C(q.s, q.v1, q.v2, q.v3)
AbstractCliffordNumber(q::Quaternion) = EvenCliffordNumber{VGA(3)}(q)

function (::Type{C})(q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return C(EvenCliffordNumber{VGA(3)}(q))
end

function convert(::Type{C}, q::Quaternion) where C<:AbstractCliffordNumber{VGA(3)} 
    return convert(C, EvenCliffordNumber{VGA(3)}(q))
end

#---Conversion to quaternions----------------------------------------------------------------------#

(::Type{H})(c::EvenCliffordNumber{VGA(3)}) where H<:Quaternion = H(Tuple(c)...)

# TODO: simplify this with functor syntax. Documenter.jl issues prevent this at the moment
"""
    Quaternion(c::AbstractCliffordNumber{VGA(3)})
    Quaternion{T}(c::AbstractCliffordNumber{VGA(3)})

Constructs a quaternion from an element of the algebra of physical space, the 3D geometric algebra
with a positive-definite signature whose even subalgebra is isomorphic to the quaternion algebra ℍ.
Any odd-grade coefficients of `c` are lost. If the type parameter `T` is supplied, the scalars of
the input are converted to type T.

If loss of odd-grade coefficients should throw an error, use `convert(Quaternion, c)` or
`convert(Quaternion{T}, c)` instead of the constructor.
"""
Quaternion(c::AbstractCliffordNumber{VGA(3)}) = Quaternion(EvenCliffordNumber{VGA(3)}(c))

function Quaternion{T}(c::AbstractCliffordNumber{VGA(3)}) where T
    return Quaternion{T}(EvenCliffordNumber{VGA(3)}(c))
end

function convert(::Type{H}, c::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
    # This fails if Quaternion would drop represented grades of c
    return H(convert(EvenCliffordNumber{VGA(3)}, c))
end

#---Promotion rules--------------------------------------------------------------------------------#

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

#---Extend slerp to some AbstractCliffordNumber types----------------------------------------------#
"""
    slerp(x::AbstractCliffordNumber{VGA(3)}, y::AbstractCliffordNumber{VGA(3)}, t::Real)
        -> EvenCliffordNumber{VGA(3)}

Performs spherical linear interpolation for rotors in the 3D vanilla geometric algebra, treating
them as if they were `Quaternion` instances.
"""
function slerp(a::AbstractCliffordNumber{VGA(3)}, b::AbstractCliffordNumber{VGA(3)}, t::Real)
    S = promote_type(scalar_type(a), scalar_type(b), typeof(t))
    (qa, qb) =  convert.(Quaternion{S}, (a, b))
    return EvenCliffordNumber{VGA(3)}(slerp(qa, qb, convert(S, t)))
end

"""
    slerp(x::Quaternion, y::AbstractCliffordNumber{VGA(3)}, t::Real) -> EvenCliffordNumber{VGA(3)}
    slerp(x::AbstractCliffordNumber{VGA(3)}, y::Quaternion, t::Real) -> EvenCliffordNumber{VGA(3)}

Performs spherical linear interpolation for rotors in the 3D vanilla geometric algebra, treating
them as if they were `Quaternion` instances. The resulting output gains Clifford number semantics,
so it is of type `EvenCliffordNumber{VGA(3)}` instead of `Quaternion`.
"""
function slerp(a::AbstractCliffordNumber{VGA(3)}, q::Quaternion, t::Real)
    S = promote_type(scalar_type(a), scalar_type(q), typeof(t))
    (qa, qb) = convert.(Quaternion{S}, (a, q))
    return EvenCliffordNumber{VGA(3)}(slerp(qa, qb, convert(S, t)))
end

function slerp(q::Quaternion, b::AbstractCliffordNumber{VGA(3)}, t::Real)
    S = promote_type(scalar_type(q), scalar_type(b), typeof(t))
    (qa, qb) = convert.(Quaternion{S}, (q, b))
    return EvenCliffordNumber{VGA(3)}(slerp(qa, qb, convert(S, t)))
end

# TODO: spherical linear interpolation for arbitrary dimensions

end
