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

isscalar(x::OddCliffordNumber) = iszero(x)
isscalar(x::KVector) = iszero(x)
isscalar(::KVector{0}) = true
isscalar(::AbstractCliffordNumber{QuadraticForm{0,0,0}}) = true
isscalar(::BaseNumber) = true

scalar(x::AbstractCliffordNumber) = x[scalar_index(x)]

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

#---Equality---------------------------------------------------------------------------------------#

function Base.:(==)(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return all(x[i] == y[i] for i in BitIndices(promote_type(typeof(x), typeof(y))))
end

# Define equality with scalar values in terms of scalar operations above
Base.:(==)(x::AbstractCliffordNumber, y::BaseNumber) = isscalar(x) && (scalar(x) == y)
Base.:(==)(x::BaseNumber, y::AbstractCliffordNumber) = isscalar(y) && (x == scalar(y))

#---Grade selection--------------------------------------------------------------------------------#
"""
    select_grade(x::CliffordNumber, g::Integer)

Returns a multivector similar to `x` where all elements not of grade `g` are equal to zero.
"""
select_grade(x::CliffordNumber, g::Integer) = typeof(x)(i -> x[i] * (count_ones(i) == g))

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
    coeff = (@inbounds x[a]) * (@inbounds y[b]) * sign_of_mult(a,b) * condition
    return C(ntuple(i -> coeff * (@inbounds(BitIndices(C)[i]) == abs(a*b)), Val(length(C))))
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
    CliffordNumbers.product_kernel(
        ::Type{T},
        x::AbstractCliffordNumber{Q},
        y::AbstractCliffordNumber{Q},
        [f = ((a,b) -> true)]
    )

Sums the products of each pair of nonzero basis blades of `x` and `y`. This can be used to to
implement various products by supplying a function `f` which acts on the indices of `x` and `y` to
return a `Bool`, and the product of the basis blades is excluded if it evaluates to `true`.
"""
@inline function product_kernel(
    ::Type{T},
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    f = ((a,b) -> true)
) where {Q,T<:AbstractCliffordNumber{Q}}
    # Tested against sum(); this is slightly better
    result = zero(T)
    for a in BitIndices(x), b in BitIndices(y)
        result += elementwise_product(T, x, y, a, b, f(a,b)::Bool)
    end
    return result
end

"""
    *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q})

Calculates the geometric product of `x` and `y`, returning the smallest type which is able to
represent all nonzero basis blades of the result.
"""
@inline function *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return product_kernel(geometric_product_type(typeof(x), typeof(y)), x, y)
end

@inline (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q}) where Q = x * y

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

function contraction_type(C1::Type{<:KVector{K1,Q}}, C2::Type{<:KVector{K2,Q}}) where {K1,K2,Q}
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
@inline function contraction(
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    ::Val{B}
) where {Q,B}
    @assert B isa Bool string(
        "Final argument must be Val(true) for left contraction or Val(false) for right contraction."
    )
    T = contraction_type(typeof(x), typeof(y))
    f = (a,b) -> (grade(a) - grade(b)) * Int8(-1)^B == grade(a*b)
    return product_kernel(T, x, y, f)
end

"""
    left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨼(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `KVector{n-m,Q}(A*B)`.
"""
function left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return contraction(x, y, Val(true))
end

"""
    right_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨽(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the right contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `KVector{m-n,Q}(A*B)`.
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
    T = contraction_type(typeof(x), typeof(y))
    f = (a,b) -> abs(grade(a) - grade(b)) == grade(a*b)
    return product_kernel(T, x, y, f)
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
@inline function wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    T = wedge_product_type(typeof(x), typeof(y))
    return product_kernel(T, x, y, has_wedge)
end

@inline wedge(x::Real, y::AbstractCliffordNumber) = x * y
@inline wedge(x::AbstractCliffordNumber, y::Real) = x * y

@inline wedge(x::BaseNumber, y::BaseNumber) = x * y

const ∧ = wedge

#---Commutator and anticommutator products---------------------------------------------------------#
"""
    commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the commutator (or antisymmetric) product, equal to `1//2 * (x*y - y*x)`,
or equivalently, `1//2 * (x*y - reverse(x*y))`.

Note that the commutator product, in general, is *not* equal to the wedge product, which may be
invoked with the `wedge` function or the `∧` operator.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    z = x*y
    return 1//2 * (z - reverse(z))
end

const × = commutator

"""
    anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the anticommutator (or symmetric) product, equal to `1//2 * (x*y + y*x)`, or
equivalently, `1//2 * (x*y + reverse(x*y))`.

Note that the dot product, in general, is *not* equal to the anticommutator product, which may be
invoked with `dot`. In some cases, the preferred operators might be the left and right contractions,
which use infix operators `⨼` and `⨽` respectively.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    z = x*y
    return 1//2 * (z + reverse(z))
end

const ⨰ = anticommutator

#---Duals------------------------------------------------------------------------------------------#
# TODO: make this work for things that aren't CliffordNumber
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
    isdegenerate(Q) && error("Cannot calculate the dual in a degenerate metric.")
    return CliffordNumber{Q}(ntuple(i -> x[dual(BitIndices(Q)[i])], Val(length(x))))
end

"""
    undual(x::CliffordNumber) -> CliffordNumber

Calculates the undual of `x`, which is equal to the left contraction of `x` with the pseudoscalar.
This function can be used to reverse the behavior of `dual()`.
"""
function undual(x::CliffordNumber{Q}) where Q
    isdegenerate(Q) && error("Cannot calculate the ;undual in a degenerate metric.")
    return CliffordNumber{Q}(ntuple(i -> x[undual(BitIndices(Q)[i])], Val(length(x))))
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
import Base: exp, ^

"""
    CliffordNumbers.exponential_type(::Type{<:AbstractCliffordNumber})
    CliffordNumbers.exponential_type(x::AbstractCliffordNumber)

Returns the type expected when exponentiating a Clifford number. This is an `EvenCliffordNumber` if
the nonzero grades of the input are even, a `CliffordNumber` otherwise.
"""
@generated function exponential_type(::Type{C}) where {Q,C<:AbstractCliffordNumber{Q}}
    T = typeof(exp(zero(numeric_type(C))))
    if all(iseven, nonzero_grades(C))
        return :(EvenCliffordNumber{Q,$T,div(elements(Q), 2)})
    else
        return :(CliffordNumber{Q,$T,elements(Q)})
    end
end

exponential_type(x::AbstractCliffordNumber) = exponential_type(typeof(x))

# KVector{K} promotes incorrectly for odd k due to broken type inference, see this issue:
# https://github.com/JuliaLang/julia/issues/53504
^(k::KVector, n::Integer) = convert(exponential_type(k), k)^n

"""
    CliffordNumbers.exp_taylor(x::AbstractCliffordNumber, order = 12)

Calculates the exponential of `x` using a Taylor expansion up to the specified order. In most cases,
12 is as sufficient number.
"""
function exp_taylor(x::AbstractCliffordNumber, order = 12)
    s = nextpow(2, abs(x))
    result = zero(exponential_type(x))
    for n in 0:order
        result += (x/s)^n / factorial(n)
    end
    return result^s
end

"""
    exp(x::AbstractCliffordNumber{Q})

Returns the natural exponential of a Clifford number.

For special cases where m squares to a scalar, the following shortcuts can be used to calculate
`exp(x)`:
  * When x^2 < 0: `exp(x) === cos(abs(x)) + x * sin(abs(x)) / abs(x)`
  * When x^2 > 0: `exp(x) === cosh(abs(x)) + x * sinh(abs(x)) / abs(x)`
  * When x^2 === 0: `exp(x) == 1 + x`

See also: [`exppi`](@ref), [`exptau`](@ref).
"""
function exp(x::AbstractCliffordNumber)
    T = exponential_type(x)
    sq = x*x
    if isscalar(sq)
        Tx = convert(T, x)
        scalar(sq) < 0 && return cos(abs(x)) + Tx * sin(abs(x)) / abs(x)
        scalar(sq) > 0 && return cosh(abs(x)) + Tx * sinh(abs(x)) / abs(x)
        return 1 + Tx
    end
    return exp_taylor(x)
end

"""
    exppi(x::AbstractCliffordNumber)

Returns the natural exponential of `π * x` with greater accuracy than `exp(π * x)` in the case where
`x^2` is a negative scalar, especially for large values of `abs(x)`.

See also: [`exp`](@ref), [`exptau`](@ref).
"""
function exppi(x::AbstractCliffordNumber)
    T = exponential_type(x)
    sq = x*x
    if isscalar(sq)
        Tx = convert(T, x)
        scalar(sq) < 0 && return cospi(abs(x)) + Tx * sinpi(abs(x)) / abs(x)
        # TODO: is this really more accurate?
        scalar(sq) > 0 && return cosh(π*abs(x)) + Tx * sinh(π*abs(x)) / abs(x)
        return 1 + x
    end
    # TODO: compare this to the regular Taylor expansion result
    return exp_taylor(π*x)
end

"""
    exptau(x::AbstractCliffordNumber)

Returns the natural exponential of `2π * x` with greater accuracy than `exp(2π * x)` in the case
where `x^2` is a negative scalar, especially for large values of `abs(x)`.

See also: [`exp`](@ref), [`exppi`](@ref).
"""
exptau(x::AbstractCliffordNumber) = exppi(2*x)
