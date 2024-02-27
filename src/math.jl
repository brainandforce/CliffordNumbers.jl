#---Equality---------------------------------------------------------------------------------------#

function Base.:(==)(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return all(x[i] == y[i] for i in BitIndices(Q))
end

#---Scalars and pseudoscalars----------------------------------------------------------------------#
"""
    isscalar(x::AbstractCliffordNumber)

Determines whether the Clifford number `x` is a scalar, meaning that all of its blades of nonzero
grade are zero.
"""
function isscalar(x::AbstractCliffordNumber)
    inds = Iterators.filter(i -> !iszero(grade(i)), BitIndices(x))
    return all(iszero, x[i] for i in inds)
end

isscalar(x::OddCliffordNumber) = iszero(x)
isscalar(x::KVector) = iszero(x)
isscalar(::KVector{0}) = true
isscalar(::AbstractCliffordNumber{QuadraticForm{0,0,0}}) = true
isscalar(::BaseNumber) = true

"""
    ispseudoscalar(m::AbstractCliffordNumber)

Determines whether the Clifford number `x` is a pseudoscalar, meaning that all of its blades with
grades below the dimension of the space are zero.
"""
function ispseudoscalar(x::AbstractCliffordNumber)
    inds = Iterators.filter(i -> grade(i) != dimension(QuadraticForm(x)), BitIndices(x))
    return all(iszero, x[i] for i in inds)
end

ispseudoscalar(x::KVector{K,Q}) where {K,Q} = (iszero(x) || K == dimension(Q))

#---Grade selection--------------------------------------------------------------------------------#
"""
    select_grade(x::CliffordNumber, g::Integer)

Returns a multivector similar to `x` where all elements not of grade `g` are equal to zero.
"""
select_grade(x::CliffordNumber, g::Integer) = typeof(x)(i -> x[i] * (count_ones(i) == g))

scalar(x::AbstractCliffordNumber{Q}) where Q = x[BitIndex{Q}()]

"""
    real(x::CliffordNumber{Q,T<:Real}) = T

Return the real (scalar) portion of a real Clifford number. 
"""
Base.real(x::AbstractCliffordNumber{Q,<:Real}) where Q = x[BitIndex{Q}()]

#---Sign changing operations-----------------------------------------------------------------------#

for f in (:reverse, :grade_involution, :conj)
    @eval begin
        function $f(x::T) where T<:AbstractCliffordNumber
            return T(ntuple(i -> x[$f(BitIndices(x)[i])], Val(length(T))))
        end
    end
end

Base.:~(x::AbstractCliffordNumber) = reverse(x)

#---Addition, negation, and subtraction------------------------------------------------------------#
import Base: +, -

+(x::T, y::T) where T<:AbstractCliffordNumber = T(Tuple(x) .+ Tuple(y))

function +(x::AbstractCliffordNumber, y::BaseNumber)
    T = promote_type(typeof(x), typeof(y))
    b = BitIndices(T)
    data = ntuple(i -> x[b[i]] + (iszero(grade(b[i])) * y), Val(length(T)))
    return T(data)
end

+(x::BaseNumber, y::CliffordNumber) = y + x

function -(x::AbstractCliffordNumber)
    data = (-).(Tuple(x))
    return similar_type(x, eltype(data))(data)
end

-(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = x + (-y)

function -(x::AbstractCliffordNumber, y::BaseNumber)
    T = promote_type(typeof(x), typeof(y))
    b = BitIndices(T)
    data = ntuple(i -> x[b[i]] - (iszero(grade(b[i])) * y), Val(length(T)))
    return T(data)
end

function -(x::BaseNumber, y::AbstractCliffordNumber)
    T = promote_type(typeof(x), typeof(y))
    b = BitIndices(T)
    data = ntuple(i -> (iszero(grade(b[i])) * x) - y[b[i]], Val(length(T)))
    return T(data)
end

# TODO: is it more efficient to define some more specific methods for some types?

#---Scalar multiplication--------------------------------------------------------------------------#
import Base.:*
import Base.:/
import Base.://

# These all have the same structure
for op in (:*, :/, ://)
    @eval begin
        function $op(x::AbstractCliffordNumber, y::BaseNumber)
            data = $op.(x.data, y)
            return similar_type(typeof(x), eltype(data))(data)
        end
        function $op(x::BaseNumber, y::AbstractCliffordNumber)
            data = $op.(x, y.data)
            return similar_type(typeof(y), eltype(data))(data)
        end
    end
end

(x::AbstractCliffordNumber)(y::BaseNumber) = x * y

#---Geometric product-----------------------------------------------------------------------------#
import Base.:*

"""
    CliffordNumbers.elementwise_product(
        [::Type{C},]
        x::AbstractCliffordNumber{Q},
        y::AbstractCliffordNumber{Q},
        a::BitIndex{Q},
        b::BitIndex{Q},
        [condition = true]
    ) where {Q,T<:AbstractCliffordNumber{Q}} -> C

Calculates the geometric product between the element of `x` indexed by `a` and the element of `y`
indexed by `b`. The result is returned as type `C`, but this can be inferred automatically if not
provided.

An optional boolean condition can be provided, which simplifies the implementation of certain
products derived from the geometric product.
"""
@inline function elementwise_product(
    ::Type{C},
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    a::BitIndex{Q},
    b::BitIndex{Q},
    condition::Bool = true
) where {Q,C<:AbstractCliffordNumber{Q}}
    inds = BitIndices(C)
    coeff = x[a] * y[b] * sign_of_mult(a,b) * condition
    return C(ntuple(i -> coeff * (inds[i] == abs(a*b)), Val(length(C))))
end

@inline function elementwise_product(
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    a::BitIndex{Q},
    b::BitIndex{Q},
    condition::Bool = true
) where Q
    return elementwise_product(promote_type(typeof(x), typeof(y)), x, y, a, b, condition)
end

@generated function geometric_product_type(
    ::Type{C1},
    ::Type{C2}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    c1_odd = all(isodd, nonzero_grades(C1))
    c2_odd = all(isodd, nonzero_grades(C2))
    c1_even = all(iseven, nonzero_grades(C1))
    c2_even = all(iseven, nonzero_grades(C2))
    P = (c1_odd && c2_even) || (c1_even && c2_odd)
    T = promote_numeric_type(C1,C2)
    if (!c1_odd && !c1_even) || (!c2_odd && !c2_even)
        return :(CliffordNumber{Q,$T,elements(Q)})
    else
        return :(Z2CliffordNumber{$P,Q,$T,div(elements(Q), 2)})
    end
end

"""
    *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q})

Calculates the geometric product of `x` and `y`, returning the smallest type which is able to
represent all nonzero basis blades of the result.
"""
function *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    T = geometric_product_type(typeof(x), typeof(y))
    return sum(elementwise_product(T, x, y, a, b) for a in eachindex(x), b in eachindex(y))
end

(x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q}) where Q = x * y

#---Scalar products--------------------------------------------------------------------------------#
"""
    scalar_product(x::CliffordNumber{Q,T1}, y::CliffordNumber{Q,T2}) -> promote_type(T1,T2)

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to an `AbstractCliffordNumber`.
"""
function scalar_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    # Only iterate through a minimal set of indices, known to be nonzero
    inds = eachindex(promote_type(typeof(x), typeof(y)))
    return sum(x[i] * y[i] * sign_of_mult(i) for i in inds)
end

"""
    abs2(x::AbstractCliffordNumber{Q,T}) -> T

Calculates the squared norm of `x`, equal to `scalar_product(x, ~x)`.
"""
Base.abs2(x::AbstractCliffordNumber) = scalar_product(x, ~x)

"""
    abs2(x::CliffordNumber{Q,T}) -> Union{Real,Complex}

Calculates the norm of `x`, equal to `sqrt(scalar_product(x, ~x))`.
"""
Base.abs(x::AbstractCliffordNumber) = sqrt(abs2(x))

"""
    normalize(x::AbstractCliffordNumber{Q}) -> AbstractCliffordNumber{Q}

Normalizes `x` so that its magnitude (as calculated by `abs2(x)`) is 1.
"""
normalize(x::AbstractCliffordNumber) = x / abs(x)

#---Contractions-----------------------------------------------------------------------------------#
#= Old implementations, left here for reference
function left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(numeric_type(x), numeric_type(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        coeff = x[a] * y[b] * sign_of_mult(a,b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (grade(b) - grade(a) == grade(a*b))
        result += CliffordNumber{Q,T}(i -> coeff * (i == (a*b).blade))
    end
    return result
end

function right_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    T = promote_type(numeric_type(x), numeric_type(y))
    result = zero(CliffordNumber{Q,T})
    for a in eachindex(x), b in eachindex(y)
        coeff = x[a] * y[b] * sign_of_mult(a,b)
        # Set to zero if the grade difference of b and a is not equal the grade of the new index
        coeff *= (grade(a) - grade(b) == grade(a*b))
        result += CliffordNumber{Q,T}(i -> coeff * (i == (a*b).blade))
    end
    return result
end
=#

function contraction_type(::Type{<:KVector{K1,Q}}, ::Type{<:KVector{K2,Q}}) where {K1,K2,Q}
    K = abs(K1 - K2) # this works in all cases, if K2 > K1 then the values are just zero
    return KVector{K,Q,promote_numeric_type(C1, C2),binomial(dimension(Q), K)}
end

function contraction_type(
    ::Type{C1},
    ::Type{C2}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    return geometric_product_type(C1, C2)
end

"""
    CliffordNumbers.contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}, ::Val)

Generic implementation of left and right contractions as well as dot products. The left contraction
is calculated if the final argument is `Val(true)`; the right contraction is calcluated if the final
argument is `Val(false)`, and the dot product is calculated for any other `Val`.

In general, code should never refer to this method directly; use `left_contraction`,
`right_contraction`, or `dot` if needed.
"""
function contraction(
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    ::Val{B}
) where {Q,B}
    T = contraction_type(typeof(x), typeof(y))
    itr = Iterators.filter(Iterators.product(eachindex(x), eachindex(y))) do t
        gdiff = grade(a) - grade(b)
        return gdiff * ifelse(B isa Bool, sign(gdiff), (-1)^B) == grade(a*b)
    end
    return sum(elementwise_product(T, x, y, a, b) for (a,b) in itr)
end

"""
    left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨼(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `grade_select(A*B, n-m)`.
"""
function left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return contraction(x, y, Val(true))
end

"""
    right_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨽(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the right contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `grade_select(A*B, m-n)`.
"""
function right_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return contraction(x, y, Val(false))
end

"""
    dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the dot product of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction when `n >= m`.
"""
function dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q 
    return contraction(x, y, Val(nothing))
end

const ⨼ = left_contraction
const ⨽ = right_contraction

"""
    hestenes_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Returns the Hestenes product: this is equal to the dot product given by `dot(x, y)` but is equal to
to zero when either `x` or `y` is a scalar.

This product is generally understood to lack utility; left and right contractions are preferred over
this product in almost every case. It is implemented for the sake of completeness.
"""
function hestenes_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return ifelse(isscalar(x) || isscalar(y), zero(promote_type(typeof(x), typeof(y))), dot(x,y))
end

#---Wedge (outer) product--------------------------------------------------------------------------#

function wedge_product_type(C1::Type{<:KVector{K1,Q}}, C2::Type{<:KVector{K2,Q}}) where {K1,K2,Q}
    return KVector{K1+K2,Q,promote_numeric_type(C1, C2),binomial(dimension(Q), K1+K2)}
end

function wedge_product_type(
    ::Type{C1},
    ::Type{C2}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    return geometric_product_type(C1, C2)
end

"""
    wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ∧(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the wedge (outer) product of two Clifford numbers `x` and `y` with quadratic form `Q`.

Note that the wedge product, in general, is *not* equal to the commutator product (or antisymmetric
product), which may be invoked with the `commutator` function or the `×` operator.
"""
function wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    T = wedge_product_type(typeof(x), typeof(y))
    # Only iterate through elements that have nonzero wedge products
    itr = Iterators.filter(t -> has_wedge(t...), Iterators.product(eachindex(x), eachindex(y)))
    return sum(elementwise_product(T, x, y, a, b) for (a,b) in itr)
end

wedge(x::Real, y::AbstractCliffordNumber) = x * y
wedge(x::AbstractCliffordNumber, y::Real) = x * y

wedge(x::BaseNumber, y::BaseNumber) = x * y

const ∧ = wedge

#---Commutator and anticommutator products---------------------------------------------------------#
"""
    commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the commutator product, equal to `1//2 * (x*y - y*x)`, or equivalently, 
`1//2 * (x*y - reverse(x*y))`.

Note that the wedge product, in general, is *not* equal to the wedge product, which may be invoked
with the `wedge` function or the `∧` operator.
"""
function commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    z = x*y
    return (x*y - reverse(x*y)) // 2
end

const × = commutator

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

#---Sandwich product-------------------------------------------------------------------------------#
"""
    sandwich(x::CliffordNumber{Q}, y::CliffordNumber{Q})

Calculates the sandwich product of `x` with `y`: `~y * x * y`, but with corrections for numerical
stability. 
"""
function sandwich(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    meat = normalize(x)
    bread = normalize(y)
    return abs(x) * ~bread * meat * bread
end

sandwich(x::BaseNumber, ::CliffordNumber) = x

#---Exponentials-----------------------------------------------------------------------------------#
"""
    CliffordNumbers.intlog2(x::Real) -> Int

Efficiently approximates the integer portion of the base-2 logarithm of `abs(x)` for calculating the
multivector exponential.

Note that this function returns `0` when `x == 0`.
"""
# log2(sqrt(2)) ≈ 0.5, so use that as a natural cutoff
intlog2(x::AbstractFloat) = last(Base.Math.frexp(sqrt(2) * x)) - !iszero(x)
intlog2(x::Real) = intlog2(Float64(x))
# It's likely not worth specializing on integers when considering the sqrt(2) factor
# intlog2(x::Integer) = 8*sizeof(x) - leading_zeros(abs(x) >>> 1)

"""
    exp(x::CliffordNumber{Q}) -> CliffordNumber{Q,<:AbstractFloat}

Returns the natural exponential of a Clifford number.

For special cases where m squares to a scalar, the following shortcuts can be used to calculate
`exp(x)`:
  * When x^2 < 0: `exp(x) === cos(abs(x)) + x * sin(abs(x)) / abs(x)`
  * When x^2 > 0: `exp(x) === cosh(abs(x)) + x * sinh(abs(x)) / abs(x)`
  * When x^2 === 0: `exp(x) == 1 + x`

See also: [`exppi`](@ref), [`exptau`](@ref).
"""
function Base.exp(x::CliffordNumber)
    sq = x^2
    if isscalar(sq)
        scalar(sq) < 0 && return cos(abs(x)) + x * sin(abs(x)) / abs(x)
        scalar(sq) > 0 && return cosh(abs(x)) + x * sinh(abs(x)) / abs(x)
        return 1 + x
    end
    # General case: Taylor expansion
    # Divide x by s which is chosen to keep its norm reasonably close to unity
    # But it's also kept as a power of 2 to make the power calculation more efficient
    s = 2^intlog2(abs(x))
    return sum((x/s)^n / factorial(n) for n in 0:12)^s
end

"""
    exppi(x::CliffordNumber)

Returns the natural exponential of `π * x` with greater accuracy than `exp(π * x)` in the case where
`x^2` is a negative scalar.

See also: [`exp`](@ref), [`exptau`](@ref).
"""
function exppi(x::CliffordNumber)
    sq = x^2
    if isscalar(sq)
        scalar(sq) < 0 && return cospi(abs(x)) + x * sinpi(abs(x)) / abs(x)
        # TODO: is this really more accurate?
        scalar(sq) > 0 && return cosh(π*abs(x)) + x * sinh(π*abs(x)) / abs(x)
        return 1 + x
    end
    # General case: Taylor expansion
    # Divide x by s which is chosen to keep its norm reasonably close to unity
    # But it's also kept as a power of 2 to make the power calculation more efficient
    s = 2^intlog2(π*abs(x))
    return sum((x/s)^n * pi^n / factorial(n) for n in 0:20)^s
end

"""
    exptau(x::CliffordNumber)

Returns the natural exponential of `2π * x` with greater accuracy than `exp(2π * x)` in the case
where `x^2` is a negative scalar.

See also: [`exp`](@ref), [`exppi`](@ref).
"""
exptau(x::CliffordNumber) = exppi(2*x)
