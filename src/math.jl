#---Grade selection--------------------------------------------------------------------------------#
"""
    select_grade(x::CliffordNumber, g::Integer)

Returns a multivector similar to `x` where all elements not of grade `g` are equal to zero.
"""
select_grade(x::CliffordNumber, g::Integer) = typeof(x)(i -> x[i] * (hamming_weight(i) == g))

scalar(x::CliffordNumber) = first(x.data)

"""
    real(x::CliffordNumber{Q,T<:Real}) = T

Return the real (scalar) portion of a real Clifford number. 
"""
Base.real(x::CliffordNumber{<:QuadraticForm,<:Real}) = first(x.data)

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

"""
    CliffordAlgebra.elementwise_product(
        x::CliffordNumber{Q},
        y::CliffordNumber{Q},
        a::BitIndex{Q},
        b::BitIndex{Q}
    )

Calculates the geometric product between the element of `x` indexed by `a` and the element of `y`
indexed by `b`.
"""
@inline function elementwise_product(
    x::CliffordNumber{Q},
    y::CliffordNumber{Q},
    a::BitIndex{Q},
    b::BitIndex{Q}
) where {Q}
    return CliffordNumber{Q}(i -> x[a] * y[b] * sign_of_mult(a,b) * (i == (a*b).blade))
end

"""
    *(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the geometric product between multivectors/Clifford numbers `x` and `y` which share the
quadratic form `Q`.
"""
function *(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        result += elementwise_product(x, y, a, b)
    end
    return result
end

#---Scalar products--------------------------------------------------------------------------------#

"""
    scalar_product(x::CliffordNumber{Q,T1}, y::CliffordNumber{Q,T2}) -> promote_type(T1,T2)

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to a `CliffordNumber`.

This is equal to `grade_select(x*y, 0)` but is significantly more efficient.
"""
function scalar_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    return sum(x[i] * y[i] * sign_of_mult(i) for i in eachindex(CliffordNumber{Q}))
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

Normalizes `x` so that its magnitude (as calculated by `abs2(x)`) is 1.
"""
normalize(x::CliffordNumber) = x / abs(x)

#---Contractions-----------------------------------------------------------------------------------#
"""
    left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `grade_select(A*B, n-m)`.
"""
function left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        coeff = x[a] * y[b] * sign_of_mult(a,b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (grade(b) - grade(a) == grade(a*b))
        result += CliffordNumber{Q,T}(i -> coeff * (i == (a*b).blade))
    end
    return result
end

"""
    right_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the right contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `grade_select(A*B, m-n)`.
"""
function right_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        coeff = x[a] * y[b] * sign_of_mult(a,b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (grade(a) - grade(b) == grade(a*b))
        result += CliffordNumber{Q,T}(i -> coeff * (i == (a*b).blade))
    end
    return result
end

"""
    dot(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the dot product of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction when `n >= m`.
"""
function dot(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        coeff = x[a] * y[b] * sign_of_mult(a,b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (abs(grade(b) - grade(a)) == grade(a*b))
        result += CliffordNumber{Q,T}(i -> coeff * (i == (a*b).blade))
    end
    return result
end

"""
    hestenes_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Returns the Hestenes product: this is equal to the dot product given by `dot(x, y)` but is equal to
to zero when either `x` or `y` is a scalar.
"""
function hestenes_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    return isscalar(x) || isscalar(y) ? zero(CliffordNumber{Q,T}) : dot(x,y)
end

#---Wedge (outer) product--------------------------------------------------------------------------#
"""
    wedge(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the wedge (outer) product of two Clifford numbers `x` and `y` with quadratic form `Q`.
"""
function wedge(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(eltype(x), eltype(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        result += elementwise_product(x, y, a, b) * iszero(a.blade & b.blade)
    end
    return result
end

wedge(x::Real, y::CliffordNumber{Q}) where Q = CliffordNumber{Q}(x .* y.data)
wedge(x::CliffordNumber{Q}, y::Real) where Q = CliffordNumber{Q}(y .* x.data)

wedge(x::BaseNumber, y::BaseNumber) = x * y

const ∧ = wedge

#---Duals------------------------------------------------------------------------------------------#
"""
    dual(x::CliffordNumber) -> CliffordNumber

Calculates the dual of `x`, which is equal to the left contraction of `x` with the inverse of the
pseudoscalar. However, 

Note that the dual has some properties that depend on the dimension and quadratic form:
  * The inverse of the unit pseudoscalar depends on the dimension of the space. Therefore, the
periodicity of 
  * If the metric is degenerate, the dual is not unique.
"""
function dual(x::CliffordNumber{Q}) where Q
    isdegenerate(Q) || return CliffordNumber{Q}(i -> x[undual(BitIndices(Q)[i])])
end

"""
    undual(x::CliffordNumber) -> CliffordNumber

Calculates the undual of `x`, which is equal to the left contraction of `x` with the pseudoscalar.
This function can be used to reverse the behavior of `dual()`.
"""
function undual(x::CliffordNumber{Q}) where Q
    isdegenerate(Q) || return CliffordNumber{Q}(i -> x[undual(BitIndices(Q)[i])])
end

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
    # Special cases: x^2 is a scalar
    sq = x^2
    if isscalar(sq)
        real(sq) < 0 && return cos(abs(x)) + x * sin(abs(x)) / abs(x)
        real(sq) > 0 && return cosh(abs(x)) + x * sinh(abs(x)) / abs(x)
        return 1 + x
    end
    # General case: Taylor expansion
    # Divide m by s which is chosen to keep norm(m) reasonably close to unity
    # But it's also kept as a power of 2 to make the power calculation more efficient
    s = Int(exp2(round(Int, log2(abs(x)))))
    return sum((x/s)^n / factorial(n) for n in 0:12)^s
end
