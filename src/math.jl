#---Grade selection--------------------------------------------------------------------------------#
"""
    select_grade(x::CliffordNumber, g::Integer)

Returns a multivector similar to `x` where all elements not of grade `g` are equal to zero.
"""
select_grade(x::CliffordNumber, g::Integer) = typeof(x)(i -> x[i] * (hamming_weight(i) == g))

#---Sign changing operations-----------------------------------------------------------------------#
"""
    reverse(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}
    ~(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculate the reverse of Clifford number `x`. This effectively reverses the products that form the
basis blades, or in other words, reverses the order of the geometric product that resulted in `x`.
"""
Base.reverse(x::CliffordNumber) = typeof(x)(i -> x[i] * Int8(-1)^!iszero(hamming_weight(i) & 2))

"""
    reverse(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}
    ~(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculate the reverse of a Clifford number. This effectively reverses the products that form the
basis blades, or in other words, reverses the order of the geometric product that resulted in `x`.
"""
Base.:~(x::CliffordNumber) = reverse(x)

"""
    grade_involution(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculates the grade involution of Clifford number `x`. This effectively multiplies all of the basis
vectors of the space by -1, which makes elements of odd grade flip sign.
"""
grade_involution(x::CliffordNumber) = typeof(x)(i -> x[i] * Int8(-1)^!iseven(hamming_weight(i)))

"""
    conj(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculates the Clifford conjugate of a Clifford number `x`. This is equal to
`grade_involution(reverse(x))`.
"""
Base.conj(x::CliffordNumber) = typeof(x)(i -> x[i] * Int8(-1)^!iszero(i+1 & 2))

#---Addition---------------------------------------------------------------------------------------#
import Base.:+

+(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q = CliffordNumber{Q}(x.data .+ y.data)

function +(x::CliffordNumber{Q}, y::BaseNumber) where Q
    return CliffordNumber{Q}(ntuple(i -> x.data[i] + (isone(i) * y), Val{length(x)}()))
end

# Adding imaginary numbers to elements of real Clifford algebras (geometric algebras) should add
# the real part to the scalar and the imaginary part to the pseudoscalar
function +(x::CliffordNumber{<:QuadraticForm,<:Real}, y::Complex)
    L = length(x)
    data = ntuple(i -> x.data[i] + (isone(i) * real(y)) + ((i == L) * imag(y)), Val{L}())
    return typeof(x)(data)
end

+(x::BaseNumber, y::CliffordNumber) = y + x

#---Negation and subtraction-----------------------------------------------------------------------#
import Base.:-

-(x::CliffordNumber{Q}) where Q = CliffordNumber{Q}((-).(x.data))
-(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q = x + (-y)

# Automatically promote 

#---Scalar multiplication--------------------------------------------------------------------------#
import Base.:*
import Base.:/
import Base.://

# These all have the same structure

for op in (:*, :/, ://)
    @eval begin
        $op(x::CliffordNumber{Q}, y::BaseNumber) where Q = CliffordNumber{Q}($op.(x.data, y))
        $op(x::BaseNumber, y::CliffordNumber{Q}) where Q = CliffordNumber{Q}($op.(x, y.data))
    end
end

#---Geometric product-----------------------------------------------------------------------------#

import Base.:*

function metric_sign(::Type{QuadraticForm{P,Q,R}}, i1::Integer, i2::Integer) where {P,Q,R}
    # For positive-definite metrics, just return 1
    iszero(Q) && iszero(R) && return Int8(1)
    # Get integers with binary 1s representing dimensions which square to negative or zero values
    q = sum(2^n for n in P .+ (1:Q); init=0)
    r = sum(2^n for n in (P + Q) .+ 1:R; init=0)
    # If any of the dimensions square to zero, return zero
    !iszero(xor(i1, i2) & r) && return Int8(0)
    # Otherwise check the number of dimensions that square to -1
    return Int8(-1)^!isevil(xor(i1, i2))
end

metric_sign(i1::BitIndex{Q}, i2::BitIndex{Q}) where Q = metric_sign(Cl, i1.i, i2.i)
metric_sign(Q::Type{<:QuadraticForm}, i::Integer) = metric_sign(Q, i, i)
metric_sign(i::BitIndex{Q}) where Q = metric_sign(Q, i.i, i.i)

"""
    CliffordAlgebra.elementwise_product(
        x::CliffordNumber{Q},
        y::CliffordNumber{Q},
        a::Integer,
        b::Integer
    )

Calculates the geometric product between element `i1` of Clifford number `m1` and element `i2` of
Clifford number `m2`.
"""
@inline function elementwise_product(
    x::CliffordNumber{Q},
    y::CliffordNumber{Q},
    a::Integer,
    b::Integer
) where {Q}
    coeff = x[a] * y[b] * sign_of_mult(a, b) * metric_sign(Q, a, b)
    return CliffordNumber{Q}(i -> coeff * (i == xor(a, b)))
end

"""
    *(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the geometric product between multivectors/Clifford numbers `x` and `y` which share the
quadratic form `Q`.
"""
function *(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for a in R, b in R
        result += elementwise_product(x, y, a, b)
    end
    return result
end

#---Scalar products--------------------------------------------------------------------------------#

"""
    scalar_product(x::CliffordNumber{Q,T1}, y::CliffordNumber{Q,T2}) -> promote_type(T1,T2)

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to a `CliffordNumber`.

This is equal to `grade_select(m1*m2, 0)` but is significantly more efficient.
"""
function scalar_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    return sum(x[i] * y[i] * sign_of_mult(i) * metric_sign(Q,i) for i in 0:elements(Q)-1)
end

"""
    abs2(x::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the squared norm of `x`, equal to `scalar_product(x, ~x)`.
"""
Base.abs2(x::CliffordNumber) = scalar_product(x, ~x)

"""
    abs2(x::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the norm of `x`, equal to `sqrt(scalar_product(x, ~x))`.
"""
Base.abs(x::CliffordNumber) = sqrt(abs2(x))

"""
    normalize(x::CliffordNumber{Q}) -> CliffordNumber{Q}

Normalizes `x` so that its magnitude is 1.
"""
normalize(x::CliffordNumber) = x / abs(x)

"""
    left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `grade_select(A*B, n-m)`.
"""
function left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for a in R, b in R
        coeff = x[a] * y[b] * sign_of_mult(a, b) * metric_sign(Q, a, b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (hamming_weight(b) - hamming_weight(a) == hamming_weight(xor(a,b)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(a,b)))
    end
    return result
end

"""
    right_contraction(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the right contraction of `m1` and `m2`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `grade_select(A*B, m-n)`.
"""
function right_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for a in R, b in R
        coeff = x[a] * y[b] * sign_of_mult(a, b) * metric_sign(Q, a, b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (hamming_weight(a) - hamming_weight(b) == hamming_weight(xor(a,b)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(a,b)))
    end
    return result
end

"""
    dot(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the dot product of `m1` and `m2`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction when `n >= m`.
"""
function dot(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for a in R, b in R
        coeff = x[a] * y[b] * sign_of_mult(a, b) * metric_sign(Q, a, b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (abs(hamming_weight(b) - hamming_weight(a)) == hamming_weight(xor(a,b)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(a,b)))
    end
    return result
end

"""
    hestenes_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Returns the Hestenes product: this is equal to the dot product given by `dot(m1, m2)` but is equal
to zero when either `m1` or `m2` is a scalar.
"""
function hestenes_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    return isscalar(x) || isscalar(y) ? zero(CliffordNumber{Q,T}) : dot(x, y)
end

#---Wedge (outer) product--------------------------------------------------------------------------#

"""
    wedge(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the wedge (outer) product of two Clifford numbers `x` and `y` with quadratic form `Q`.
"""
function wedge(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for a in R, b in R
        result += elementwise_product(x, y, a, b) * iszero(a & b)
    end
    return result
end

wedge(x::Real, y::CliffordNumber{Q}) where Q = CliffordNumber{Q}(x .* y.data)
wedge(x::CliffordNumber{Q}, y::Real) where Q = CliffordNumber{Q}(y .* x.data)

wedge(x::BaseNumber, y::BaseNumber) = x * y

const ∧ = wedge

"""
    ⋆(x::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the Hodge dual of `x`, equivalent to multiplying `x` by its corresponding pseudoscalar.
"""
⋆(x::CliffordNumber) = x * pseudoscalar(x)

#---Division---------------------------------------------------------------------------------------#
"""
    versor_inverse(x::CliffordNumber)

Calculates the versor inverse of `x`, equal to `x / scalar_product(x, ~x)`, so that
`x * inv(x) == inv(x) * x == 1`.

The versor inverse is only guaranteed to be an inverse for blades and versors. Not all Clifford
numbers have a well-defined inverse, since Clifford numbers have zero divisors (for instance, in the
algebra of physical space, 1 + e₁ has a zero divisor).
"""
versor_inverse(x::CliffordNumber) = ~x / abs2(x)

#---Exponentials-----------------------------------------------------------------------------------#
"""
    exp(x::CliffordNumber{Q}) -> CliffordNumber{Q,<:AbstractFloat}

Returns the natural exponential of a Clifford number.

For special cases where m squares to a scalar, the following shortcuts can be used to calculate
`exp(x)`:
"""
function Base.exp(x::CliffordNumber)
    # Special cases: m^2 is a scalar
    if isscalar(x^2)
        iszero(x^2) && return 1 + x
        x^2 < 0 && return cos(abs(x)) + x * sin(abs(x)) / abs(x)
        x^2 > 0 && return cosh(abs(x)) + x * sinh(abs(x)) / abs(x)
    end
    # General case: Taylor expansion
    # Divide m by s which is chosen to keep norm(m) reasonably close to unity
    # But it's also kept as a power of 2 to make the power calculation more efficient
    s = Int(exp2(round(Int, log2(abs(x)))))
    return sum((x/s)^n / factorial(n) for n in 0:12)^s
end
