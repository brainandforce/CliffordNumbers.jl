#---Addition--------------------------------------------------------------------------------------#
import Base.:+

+(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q = CliffordNumber{Q}(m1.data .+ m2.data)

function +(m::CliffordNumber{Q}, n::Union{Real,Complex}) where Q
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

#---Negation and subtraction----------------------------------------------------------------------#
import Base.:-

-(m::CliffordNumber{Q}) where Q = CliffordNumber{Cl}((-).(m.data))
-(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q = m1 + (-m2)

# Automatically promote 

#---Scalar multiplication-------------------------------------------------------------------------#
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

#---Dot (inner) product---------------------------------------------------------------------------#

"""
    dot(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> Number
    m1 · m2 -> CliffordNumber{Q}

Calculates the dot (inner) product of two Clifford numbers with quadratic form `Cl`. The result is a
`Real` or `Complex` number. This can be converted back to a `CliffordNumber`.
"""
function dot(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q
    return sum(m1[i] * m2[i] * sign_of_mult(Q, i) for i in 0:elements(Q)-1)
end

#---Wedge (outer) product-------------------------------------------------------------------------#

"""
    wedge(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) -> CliffordNumber{Q}

Calculates the wedge (outer) product of two Clifford numbers with quadratic form `Q`. The result
is another `CliffordNumber{Q}`.
"""
function wedge(m1::CliffordNumber{Q}, m2::CliffordNumber{Q}) where Q

end

const ∧ = wedge

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

"""
    ⋆(m::CliffordNumber) -> CliffordNumber

Calculates the Hodge dual of `m`, equivalent to multiplying `m` by its corresponding pseudoscalar.
"""
⋆(m::CliffordNumber) = m * pseudoscalar(m)

#---Exponentials----------------------------------------------------------------------------------#
import Base: ^, exp

#=
function ^(m::CliffordNumber, x::Number)

end
=#
