#---Grade selection--------------------------------------------------------------------------------#
"""
    select_grade(m::CliffordNumber, g::Integer)

Returns a multivectors where all elements not of grade `g` are equal to zero.
"""
select_grade(m::CliffordNumber, g::Integer) = typeof(m)(i -> m[i] * (hamming_weight(i) == g))

#---Sign changing operations-----------------------------------------------------------------------#
"""
    reverse(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}
    ~(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculate the reverse of a Clifford number. This effectively reverses the products that form the
basis blades.
"""
Base.reverse(m::CliffordNumber) = typeof(m)(i -> m[i] * Int8(-1)^!iszero(hamming_weight(i) & 2))

"""
    reverse(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}
    ~(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculate the reverse of a Clifford number. This effectively reverses the products that form the
basis blades, or in other words, reverses the order of the geometric product that resulted in `m`.
"""
Base.:~(m::CliffordNumber) = reverse(m)

"""
    grade_involution(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculates the grade involution of a Clifford number. This effectively multiplies all of the basis
vectors of the space by -1, which makes elements of odd grade flip sign.
"""
grade_involution(m::CliffordNumber) = typeof(m)(i -> m[i] * Int8(-1)^!iseven(hamming_weight(i)))

"""
    conj(m::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}

Calculates the Clifford conjugate of a Clifford number. This is equal to
`grade_involution(reverse(m))`.
"""
Base.conj(m::CliffordNumber) = typeof(m)(i -> m[i] * Int8(-1)^!iszero(i+1 & 2))

#---Addition---------------------------------------------------------------------------------------#
import Base.:+

+(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q = CliffordNumber{Q}(m1.data .+ m2.data)

function +(m::CliffordNumber{Q}, n::BaseNumber) where Q
    return CliffordNumber{Q}(ntuple(i -> m.data[i] + (isone(i) * n), Val{length(m)}()))
end

# Adding imaginary numbers to elements of real Clifford algebras (geometric algebras) should add
# the real part to the scalar and the imaginary part to the pseudoscalar
function +(m::CliffordNumber{<:QuadraticForm,<:Real}, n::Complex)
    L = length(m)
    data = ntuple(i -> m.data[i] + (isone(i) * real(n)) + ((i == L) * imag(n)), Val{L}())
    return typeof(m)(data)
end

+(n::BaseNumber, m::CliffordNumber) = m + n

#---Negation and subtraction-----------------------------------------------------------------------#
import Base.:-

-(m::CliffordNumber{Q}) where Q = CliffordNumber{Q}((-).(m.data))
-(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q = m1 + (-m2)

# Automatically promote 

#---Scalar multiplication--------------------------------------------------------------------------#
import Base.:*
import Base.:/
import Base.://

# These all have the same structure

for op in (:*, :/, ://)
    @eval begin
        $op(m::CliffordNumber{Q}, n::BaseNumber) where Q = CliffordNumber{Q}($op.(m.data, n))
        $op(n::BaseNumber, m::CliffordNumber{Q}) where Q = CliffordNumber{Q}($op.(n, m.data))
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
        m1::CliffordNumber{Q},
        m2::CliffordNumber{Q},
        i1::Integer,
        i2::Integer
    )

Calculates the geometric product between element `i1` of Clifford number `m1` and element `i2` of
Clifford number `m2`.
"""
@inline function elementwise_product(
    m1::CliffordNumber{Q},
    m2::CliffordNumber{Q},
    i1::Integer,
    i2::Integer
) where {Q}
    coeff = m1[i1] * m2[i2] * sign_of_mult(i1, i2) * metric_sign(Q, i1, i2)
    return CliffordNumber{Q}(i -> coeff * (i == xor(i1, i2)))
end

"""
    *(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the geometric product between multivectors/Clifford numbers `m1` and `m2` which share the
quadratic form `Q`.
"""
function *(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for i1 in R, i2 in R
        result += elementwise_product(m1, m2, i1, i2)
    end
    return result
end

#---Scalar products--------------------------------------------------------------------------------#

"""
    scalar_product(m1::CliffordNumber{Q,T1}, m2::CliffordNumber{Q,T2}) -> promote_type(T1,T2)

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to a `CliffordNumber`.

This is equal to `grade_select(m1*m2, 0)` but is significantly more efficient.
"""
function scalar_product(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    return sum(m1[i] * m2[i] * sign_of_mult(i) * metric_sign(Q,i) for i in 0:elements(Q)-1)
end

"""
    abs2(m::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the squared norm of `m`, equal to `scalar_product(m, ~m)`.
"""
Base.abs2(m::CliffordNumber) = scalar_product(m, ~m)

"""
    abs2(m::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the norm of `m`, equal to `sqrt(scalar_product(m, ~m))`.
"""
Base.abs(m::CliffordNumber) = sqrt(abs2(m))

"""
    normalize(m::CliffordNumber{Q}) -> CliffordNumber{Q}

Normalizes `m` so that its magnitude is 1.
"""
normalize(m::CliffordNumber) = m / abs(m)

"""
    left_contraction(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the left contraction of `m1` and `m2`.

For basis blades `A` of grade `a` and `B` of grade `b`, the left contraction is zero if `b < a`,
otherwise it is `grade_select(A*B, b-a)`.
"""
function left_contraction(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for i1 in R, i2 in R
        coeff = m1[i1] * m2[i2] * sign_of_mult(i1, i2) * metric_sign(Q, i1, i2)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (hamming_weight(i2) - hamming_weight(i1) == hamming_weight(xor(i1,i2)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(i1,i2)))
    end
    return result
end

"""
    right_contraction(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the right contraction of `m1` and `m2`.

For basis blades `A` of grade `a` and `B` of grade `b`, the right contraction is zero if `a < b`,
otherwise it is `grade_select(A*B, a-b)`.
"""
function right_contraction(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for i1 in R, i2 in R
        coeff = m1[i1] * m2[i2] * sign_of_mult(i1, i2) * metric_sign(Q, i1, i2)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (hamming_weight(i1) - hamming_weight(i2) == hamming_weight(xor(i1,i2)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(i1,i2)))
    end
    return result
end

"""
    dot(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the dot product of `m1` and `m2`.

For basis blades `A` of grade `a` and `B` of grade `b`, the dot product is equal to the left
contraction when `a >= b` and is equal to the right contraction when `b >= a`.
"""
function dot(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for i1 in R, i2 in R
        coeff = m1[i1] * m2[i2] * sign_of_mult(i1, i2) * metric_sign(Q, i1, i2)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (abs(hamming_weight(i2) - hamming_weight(i1)) == hamming_weight(xor(i1,i2)))
        result += CliffordNumber{Q,T}(i -> coeff * (i == xor(i1,i2)))
    end
    return result
end

"""
    hestenes_product(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Returns the Hestenes product: this is equal to the dot product given by `dot(m1, m2)` but is equal
to zero when either `m1` or `m2` is a scalar.
"""
function hestenes_product(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    return isscalar(m1) || isscalar(m2) ? zero(CliffordNumber{Q,T}) : dot(m1, m2)
end

#---Wedge (outer) product--------------------------------------------------------------------------#

"""
    wedge(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the wedge (outer) product of two Clifford numbers with quadratic form `Q`. The result
is another `CliffordNumber{Q}`.
"""
function wedge(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Q) - 1
    result = zero(CliffordNumber{Q,T})
    for i1 in R, i2 in R
        result += elementwise_product(m1, m2, i1, i2) * iszero(i1 & i2)
    end
    return result
end

wedge(s::Real, m::CliffordNumber{Q}) where Q = CliffordNumber{Q}(s .* m.data)
wedge(m::CliffordNumber{Q}, s::Real) where Q = CliffordNumber{Q}(s .* m.data)

wedge(x::BaseNumber, y::BaseNumber) = x * y

const ∧ = wedge

"""
    ⋆(m::CliffordNumber) -> CliffordNumber

Calculates the Hodge dual of `m`, equivalent to multiplying `m` by its corresponding pseudoscalar.
"""
⋆(m::CliffordNumber) = m * pseudoscalar(m)

#---Division---------------------------------------------------------------------------------------#
"""
    versor_inverse(m::CliffordNumber)

Calculates the versor inverse of `m`, equal to `m / scalar_product(m, ~m)`, so that
`m * inv(m) == inv(m) * m == 1`.

The versor inverse may not always be well defined, particularly when the quadratic form is
degenerate or `m` has null factors.
"""
versor_inverse(m::CliffordNumber) = ~m / abs2(m)

#---Exponentials-----------------------------------------------------------------------------------#
import Base: ^, exp

#=
function ^(m::CliffordNumber, x::Number)

end
=#
