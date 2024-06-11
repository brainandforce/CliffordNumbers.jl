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

isscalar(x::Union{CliffordNumber,EvenCliffordNumber}) = all(iszero, Tuple(x)[2:end])
isscalar(x::OddCliffordNumber) = iszero(x)
isscalar(x::KVector) = iszero(x)
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
"""
function scalar_product(x::AbstractCliffordNumber{Q,T}, y::AbstractCliffordNumber{Q,T}) where {Q,T}
    # This is so much faster than anything else I can come up with
    return @inline scalar(mul(x, y, GradeFilter{:*}()))
end

"""
    abs2(x::AbstractCliffordNumber{Q,T}) -> T

Calculates the squared norm of `x`, equal to `scalar_product(x, x')`.
"""
abs2(x::AbstractCliffordNumber) = scalar_product(x, x')

"""
    abs(x::AbstractCliffordNumber{Q,T}) -> Union{Real,Complex}

Calculates the norm of `x`, equal to `sqrt(scalar_product(x, x'))`.
"""
abs(x::AbstractCliffordNumber) = hypot(Tuple(x)...)

"""
    normalize(x::AbstractCliffordNumber{Q}) -> AbstractCliffordNumber{Q}

Normalizes `x` so that its magnitude (as calculated by `abs2(x)`) is 1.
"""
normalize(x::AbstractCliffordNumber) = x / abs(x)

# Optimized versions for `KVector` scalars
abs(k::KVector{0}) = abs(scalar(k))
abs2(k::KVector{0}) = abs2(scalar(k))
normalize(x::KVector{0}) = one(typeof(x))

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
