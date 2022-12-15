const OtherNumber = Union{Real,Complex}

#---Addition--------------------------------------------------------------------------------------#
import Base.:+

+(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) where Cl = CliffordNumber{Cl}(m1.data .+ m2.data)

function +(m::CliffordNumber{Cl}, n::Union{Real,Complex}) where Cl
    return CliffordNumber{Cl}(ntuple(i -> m.data[i] + (isone(i) * n), Val{length(m)}()))
end

# Adding imaginary numbers to elements of real Clifford algebras (geometric algebras) should add
# the real part to the scalar and the imaginary part to the pseudoscalar
function +(m::CliffordNumber{<:QuadraticForm,<:Real}, n::Complex)
    L = length(m)
    data = ntuple(i -> m.data[i] + (isone(i) * real(n)) + ((i == L) * imag(n)), Val{L}())
    return typeof(m)(data)
end

+(n::OtherNumber, m::CliffordNumber) = m + n

#---Negation and subtraction----------------------------------------------------------------------#
import Base.:-

-(m::CliffordNumber{Cl}) where Cl = CliffordNumber{Cl}((-).(m.data))
-(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) where Cl = m1 + (-m2)

# Automatically promote 

#---Scalar multiplication-------------------------------------------------------------------------#
import Base.:*
import Base.:/
import Base.://

# These all have the same structure

for op in (:*, :/, ://)
    @eval begin
        $op(m::CliffordNumber{Cl}, n::OtherNumber) where Cl = CliffordNumber{Cl}($op.(m.data, n))
        $op(n::OtherNumber, m::CliffordNumber{Cl}) where Cl = CliffordNumber{Cl}($op.(n, m.data))
    end
end

#---Dot (inner) product---------------------------------------------------------------------------#

"""
    dot(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) -> Number
    m1 · m2 -> CliffordNumber{Cl}

Calculates the dot (inner) product of two Clifford numbers with quadratic form `Cl`. The result is
a `Real` or `Complex` number. This can be converted back to a `CliffordNumber`.
"""
function dot(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) where Cl
    return sum(m1[i] * m2[i] * sign_of_mult(Cl, i) for i in 0:elements(Cl)-1)
end

#---Wedge (outer) product-------------------------------------------------------------------------#

"""
    wedge(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) -> CliffordNumber{Cl}

Calculates the wedge (outer) product of two Clifford numbers with quadratic form `Cl`. The result
is another `CliffordNumber{Cl}`.
"""
function wedge(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) where Cl

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

metric_sign(i1::BitIndex{Cl}, i2::BitIndex{Cl}) where Cl = metric_sign(Cl, i1.i, i2.i)
metric_sign(Cl::Type{<:QuadraticForm}, i::Integer) = metric_sign(Cl, i, i)
metric_sign(i::BitIndex{Cl}) where Cl = metric_sign(Cl, i.i, i.i)

"""
    CliffordAlgebra.elementwise_product(
        m1::CliffordNumber{Cl},
        m2::CliffordNumber{Cl},
        i1::Integer,
        i2::Integer
    )

Calculates the geometric product between element `i1` of Clifford number `m1` and element `i2` of
Clifford number `m2`.
"""
@inline function elementwise_product(
    m1::CliffordNumber{Cl},
    m2::CliffordNumber{Cl},
    i1::Integer,
    i2::Integer
) where {Cl}
    coeff = m1[i1] * m2[i2] * sign_of_mult(i1, i2) * metric_sign(Cl, i1, i2)
    return CliffordNumber{Cl}(i -> coeff * (i == xor(i1, i2)))
end

"""
    *(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) -> CliffordNumber{Cl}

Calculates the geometric product between multivectors/Clifford numbers `m1` and `m2` which share
the quadratic form `Cl`.
"""
function *(m1::CliffordNumber{Cl}, m2::CliffordNumber{Cl}) where Cl
    T = promote_type(eltype(m1), eltype(m2))
    R = 0:elements(Cl) - 1
    result = zero(CliffordNumber{Cl,T})
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
