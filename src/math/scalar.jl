#---Scalars and pseudoscalars----------------------------------------------------------------------#
"""
    isscalar(x::AbstractCliffordNumber)

Determines whether the Clifford number `x` is a scalar, meaning that all of its blades of nonzero
grade are zero.
"""
function isscalar(x::AbstractCliffordNumber)
    inds = Iterators.filter(!isequal(scalar_index(x)), BitIndices(x))
    return all(iszero, x[i] for i in inds)
end

# TODO: how to test the fallback method which is not used by any internal types?

# The first index of these types is the scalar index
isscalar(x::Union{CliffordNumber,EvenCliffordNumber}) = all(iszero, Tuple(x)[2:end])
# These types cannot be scalars unless they are zero
isscalar(x::OddCliffordNumber) = iszero(x)
isscalar(x::KVector) = iszero(x)
# These types are definitely scalars
isscalar(::KVector{0}) = true
isscalar(::BaseNumber) = true

"""
    scalar(x::AbstractCliffordNumber{Q,T}) -> T

Returns the scalar portion of `x` as its scalar type. This is equivalent to `x[scalar_index(x)]`.

To retain Clifford number semantics, use the `KVector{0}` constructor.
"""
scalar(x::AbstractCliffordNumber) = x[scalar_index(x)]
scalar(x::KVector{0}) = only(Tuple(x))

"""
    ispseudoscalar(m::AbstractCliffordNumber)

Determines whether the Clifford number `x` is a pseudoscalar, meaning that all of its blades with
grades below the dimension of the space are zero.
"""
function ispseudoscalar(x::AbstractCliffordNumber)
    inds = Iterators.filter(!isequal(pseudoscalar_index(x)), BitIndices(x))
    return all(iszero, x[i] for i in inds)
end

ispseudoscalar(x::KVector{K,Q}) where {K,Q} = (iszero(x) || K == dimension(Q))

#---Scalar products--------------------------------------------------------------------------------#
"""
    scalar_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to an `AbstractCliffordNumber`.

The result is equal to `scalar(x * y)`, but does not calculate the coefficients associated with any
other basis blades.
"""
@inline @generated function scalar_product(x::T, y::T) where T<:AbstractCliffordNumber
    ind_signs = map(sign_of_square, BitIndices(T))
    return :(mapreduce(*, +, Tuple(x), Tuple(y), $ind_signs))
end

function scalar_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    # TODO: we could actually demote types here.
    # Indices not represented by both x and y should be skipped and add zero.
    return @inline scalar_product(promote(x, y)...)
end

function scalar_product(x::AbstractCliffordNumber, y::AbstractCliffordNumber)
    throw(AlgebraMismatch(scalar_product, (x, y)))
end

"""
    abs2(x::AbstractCliffordNumber{Q,T}) -> T

Calculates the modulus of `x` by calculating the scalar product of `x` with its reverse:
`scalar_product(x, x')`. In positive-definite metrics, this value is positive for any nonzero 
multivector and equal to `zero(T)` for any multivector equal to zero.
"""
abs2(x::AbstractCliffordNumber) = scalar_product(x, x')

"""
    abs(x::AbstractCliffordNumber{Q}) -> Union{Real,Complex}

Calculates the norm of `x`, equal to `sqrt(abs(abs2(x)))`.

The inclusion of `abs` in this expression accounts for the possibility that the algebra `Q` contains
1-blades with a negative square, which would result in `abs(x)` being imaginary. In these cases,
`abs(x)` may not be a norm, but it is used internally by `normalize(x)` to calculate a
normalization factor.
"""
abs(x::AbstractCliffordNumber) = sqrt(abs(abs2(x)))

"""
    normalize(x::AbstractCliffordNumber{Q}) -> AbstractCliffordNumber{Q}

Normalizes `x` so that its modulus (as calculated by `abs2`) is 1, 0, or -1. This procedure cannot
change the sign of the modulus.

If `abs2(x)` is zero (possible in any non-positive-definite metric), `x` is returned unchanged.
"""
normalize(x::AbstractCliffordNumber) = (n = abs(x); ifelse(iszero(n), x, x / n))

# Optimized versions for `KVector` scalars
abs(k::KVector{0}) = abs(scalar(k))
abs2(k::KVector{0}) = abs2(scalar(k))
normalize(x::KVector{0}) = ifelse(iszero(x), zero(x), one(x))

#---Scalar multiplication and division-------------------------------------------------------------#

@inline function *(x::AbstractCliffordNumber, y::BaseNumber)
    data = map(_x -> (_x * y), Tuple(x))
    return similar_type(typeof(x), eltype(data))(data)
end

# Don't assume commutative multiplication, just to be safe
@inline function *(x::BaseNumber, y::AbstractCliffordNumber)
    data = map(_y -> (x * _y), Tuple(y))
    return similar_type(typeof(y), eltype(data))(data)
end

for op in (:/, ://)
    @eval begin
        @inline function $op(x::AbstractCliffordNumber, y::BaseNumber)
            data = map(_x -> $op(_x, y), Tuple(x))
            return similar_type(typeof(x), eltype(data))(data)
        end
        @inline function $op(x::BaseNumber, y::AbstractCliffordNumber)
            data = Tuple(y') .* $op(x, abs2(y))
            return similar_type(typeof(y), eltype(data))(data)
        end
    end
end

# Needed for ambiguity resolution with //(x::Number, y::Complex)
function //(x::AbstractCliffordNumber, y::Complex)
    data = map(_x -> //(_x, y), tuple(x))
    return similar_type(typeof(x), eltype(data))(data)
end

(x::AbstractCliffordNumber)(y::BaseNumber) = x * y

"""
    muladd(x::Union{Real,Complex}, y::AbstractCliffordNumber{Q}, z::AbstractCliffordNumber{Q})
    muladd(x::AbstractCliffordNumber{Q}, y::Union{Real,Complex}, z::AbstractCliffordNumber{Q})

Multiplies a scalar with a Clifford number and adds another Clifford number, utilizing optimizations
made available with scalar `muladd`, such as `fma` if hardware support is available.
"""
function muladd(x::BaseNumber, y::T, z::T) where T<:AbstractCliffordNumber
    return T(map((_y, _z) -> muladd(x, _y, _z), Tuple(y), Tuple(z)))
end

function muladd(x::T, y::BaseNumber, z::T) where T<:AbstractCliffordNumber
    return T(map((_x, _z) -> muladd(_x, y, _z), Tuple(x), Tuple(z)))
end

function muladd(x::BaseNumber, y::AbstractCliffordNumber{Q}, z::AbstractCliffordNumber{Q}) where Q
    (yy, zz) = promote(y, z)
    return muladd(x, yy, zz)
end

function muladd(x::AbstractCliffordNumber{Q}, y::BaseNumber, z::AbstractCliffordNumber{Q}) where Q
    (xx, zz) = promote(x, z)
    return muladd(xx, y, zz)
end
