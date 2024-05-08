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

#---Equality and approximate equality--------------------------------------------------------------#

function ==(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return all(x[i] == y[i] for i in BitIndices(promote_type(typeof(x), typeof(y))))
end

# Define equality with scalar values in terms of scalar operations above
==(x::AbstractCliffordNumber, y::BaseNumber) = isscalar(x) && (scalar(x) == y)
==(x::BaseNumber, y::AbstractCliffordNumber) = isscalar(y) && (x == scalar(y))

function isapprox(
    x::AbstractCliffordNumber,
    y::AbstractCliffordNumber;
    atol::Real = 0,
    rtol::Real = Base.rtoldefault(scalar_type(x), scalar_type(y), atol),
    nans::Bool = false,
    norm = abs
)
    # This mirrors the implementation for `AbstractArray`
    d = abs(x - y)
    if isfinite(d)
        return d <= max(atol, rtol*max(norm(x), norm(y)))
    else
        (tx, ty) = Tuple.(promote(x, y))
        return all(isapprox.(tx, ty; atol, rtol, nans, norm))
    end
    #= TODO:
    If we strip away the `isfinite` if block to simplify the function, we get *worse* performance!
    This is probably due to the closure bug, but why does the if block avoid the problem?
    =#
end

#---Sign changing operations-----------------------------------------------------------------------#

for f in (:adjoint, :grade_involution, :conj)
    @eval $f(x::T) where T<:AbstractCliffordNumber = x[$f.(BitIndices(T))]
end

reverse(x::AbstractCliffordNumber) = adjoint(x)

# Faster implementations for KVector that don't require indexing
adjoint(x::T) where T<:KVector = T(x.data .* Int8(-1)^!iszero(grade(x) & 2))
grade_involution(x::T) where T<:KVector = T(x.data .* Int8(-1)^isodd(grade(x)))
conj(x::T) where T<:KVector = T(x.data .* Int8(-1)^!iszero((grade(x) + 1) & 2))

#---Addition, negation, and subtraction------------------------------------------------------------#

+(x::T, y::T) where T<:AbstractCliffordNumber = T(map(+, Tuple(x), Tuple(y)))

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

for op in (:*, :/, ://)
    @eval begin
        @inline function $op(x::AbstractCliffordNumber, y::BaseNumber)
            data = map(_x -> $op(_x, y), Tuple(x))
            return similar_type(typeof(x), eltype(data))(data)
        end
        @inline function $op(x::BaseNumber, y::AbstractCliffordNumber)
            data = map(_y -> $op(x, _y), Tuple(y))
            return similar_type(typeof(y), eltype(data))(data)
        end
    end
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

#---Geometric product-----------------------------------------------------------------------------#
"""
    *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q})

Calculates the geometric product of `x` and `y`, returning the smallest type which is able to
represent all nonzero basis blades of the result.
"""
@inline function *(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return mul(scalar_promote(x, y)...)
end

# KVector{0,Q} is just a scalar compatible with an AbstractCliffordNumber{Q}
# TODO: this may be best in an eval block with all other products
*(k::KVector{0,Q}, x::AbstractCliffordNumber{Q}) where Q = only(Tuple(k)) * x
*(x::AbstractCliffordNumber{Q}, k::KVector{0,Q}) where Q = x * only(Tuple(k))
*(k::KVector{0,Q}, l::KVector{0,Q}) where Q = KVector{0,Q}(only(Tuple(k)) * only(Tuple(l)))

@inline (x::AbstractCliffordNumber{Q})(y::AbstractCliffordNumber{Q}) where Q = x * y

#---Scalar products--------------------------------------------------------------------------------#
"""
    scalar_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the scalar product of two Clifford numbers with quadratic form `Q`. The result is a
`Real` or `Complex` number. This can be converted back to an `AbstractCliffordNumber`.
"""
function scalar_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    # Only iterate through a minimal set of indices, known to be nonzero
    T = promote_scalar_type(x, y)
    result = zero(T)
    inds = eachindex(promote_type(typeof(x), typeof(y)))
    for i in inds
        result += T(x[i] * y[i] * sign_of_square(i))
    end
    return result
end

"""
    abs2(x::AbstractCliffordNumber{Q,T}) -> T

Calculates the squared norm of `x`, equal to `scalar_product(x, x')`.
"""
function abs2(x::AbstractCliffordNumber)
    result = zero(scalar_type(x))
    for i in eachindex(Tuple(x))
        result += Tuple(x)[i] * Tuple(x)[i]
    end
    return result
end

"""
    abs2(x::CliffordNumber{Q,T}) -> Union{Real,Complex}

Calculates the norm of `x`, equal to `sqrt(scalar_product(x, x'))`.
"""
abs(x::AbstractCliffordNumber) = hypot(Tuple(x)...)

"""
    normalize(x::AbstractCliffordNumber{Q}) -> AbstractCliffordNumber{Q}

Normalizes `x` so that its magnitude (as calculated by `abs2(x)`) is 1.
"""
normalize(x::AbstractCliffordNumber) = x / abs(x)

#---Contractions-----------------------------------------------------------------------------------#
"""
    left_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨼(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the left contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the left contraction is zero if `n < m`,
otherwise it is `KVector{n-m,Q}(A*B)`.
"""
function ⨼(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return mul(scalar_promote(x, y)..., GradeFilter{:⨼}())
end

"""
    right_contraction(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ⨽(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the right contraction of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the right contraction is zero if `m < n`,
otherwise it is `KVector{m-n,Q}(A*B)`.
"""
function ⨽(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return mul(scalar_promote(x, y)..., GradeFilter{:⨽}())
end

"""
    dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the dot product of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction when `n >= m`.
"""
function dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q 
    return mul(scalar_promote(x, y)..., GradeFilter{:dot}())
end

const left_contraction = ⨼
const right_contraction = ⨽

"""
    hestenes_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Returns the Hestenes product: this is equal to the dot product given by `dot(x, y)` but is equal to
to zero when either `x` or `y` is a scalar.

This product is generally understood to lack utility; left and right contractions are preferred over
this product in almost every case. It is implemented for the sake of completeness.
"""
function hestenes_product(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return dot(x, y) * !(isscalar(x) || isscalar(y))
end

#---Wedge (outer) product--------------------------------------------------------------------------#
"""
    ∧(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    wedge(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the wedge (outer) product of two Clifford numbers `x` and `y` with quadratic form `Q`.

Note that the wedge product, in general, is *not* equal to the commutator product (or antisymmetric
product), which may be invoked with the `commutator` function or the `×` operator.
"""
function ∧(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return mul(scalar_promote(x, y)..., GradeFilter{:∧}())
end

∧(x::BaseNumber, y::BaseNumber) = x * y
∧(x::BaseNumber, y::AbstractCliffordNumber) = x * y
∧(x::AbstractCliffordNumber, y::BaseNumber) = y * x

const wedge = ∧

#---Commutator and anticommutator products---------------------------------------------------------#
"""
    commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the commutator (or antisymmetric) product, equal to `1//2 * (x*y - y*x)`.

Note that the commutator product, in general, is *not* equal to the wedge product, which may be
invoked with the `wedge` function or the `∧` operator.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y - y*x)
end

const × = commutator

"""
    anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the anticommutator (or symmetric) product, equal to `1//2 * (x*y + y*x)`.

Note that the dot product, in general, is *not* equal to the anticommutator product, which may be
invoked with `dot`. In some cases, the preferred operators might be the left and right contractions,
which use infix operators `⨼` and `⨽` respectively.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y + y*x)
end

const ⨰ = anticommutator

#---Duals------------------------------------------------------------------------------------------#
# TODO: make this work for things that aren't CliffordNumber
"""
    dual(x::CliffordNumber) -> CliffordNumber

Calculates the dual of `x`, which is equal to the left contraction of `x` with the inverse of the
pseudoscalar.

Note that the dual has some properties that depend on the dimension and quadratic form:
  * The inverse of the unit pseudoscalar is equal to its reverse, meaning that the sign may be
positive or negative depending on the total number of dimensions in the space.
  * If the metric is degenerate, the dual is not unique.
"""
function dual(x::CliffordNumber{Q}) where Q
    is_degenerate(Q) && error("Cannot calculate the dual in a degenerate metric.")
    return CliffordNumber{Q}(ntuple(i -> x[dual(BitIndices(Q)[i])], Val(length(x))))
end

"""
    undual(x::CliffordNumber) -> CliffordNumber

Calculates the undual of `x`, which is equal to the left contraction of `x` with the pseudoscalar.
This function can be used to reverse the behavior of `dual()`.
"""
function undual(x::CliffordNumber{Q}) where Q
    is_degenerate(Q) && error("Cannot calculate the ;undual in a degenerate metric.")
    return CliffordNumber{Q}(ntuple(i -> x[undual(BitIndices(Q)[i])], Val(length(x))))
end

#---Division---------------------------------------------------------------------------------------#
"""
    versor_inverse(x::CliffordNumber)

Calculates the versor inverse of `x`, equal to `x / scalar_product(x, x')`, so that
`x * inv(x) == inv(x) * x == 1`.

The versor inverse is only guaranteed to be an inverse for blades and versors. Not all Clifford
numbers have a well-defined inverse, since Clifford numbers have zero divisors (for instance, in the
algebra of physical space, 1 + e₁ has a zero divisor).
"""
versor_inverse(x::CliffordNumber) = x' / abs2(x)

#---Sandwich product-------------------------------------------------------------------------------#
"""
    sandwich(x::CliffordNumber{Q}, y::CliffordNumber{Q})

Calculates the sandwich product of `x` with `y`: `y' * x * y`, but with corrections for numerical
stability. 
"""
function sandwich(x::CliffordNumber{Q}, y::CliffordNumber{Q}) where Q
    meat = normalize(x)
    bread = normalize(y)
    return abs(x) * bread' * meat * bread
end

sandwich(x::BaseNumber, ::CliffordNumber) = x

#---Exponentials-----------------------------------------------------------------------------------#
"""
    CliffordNumbers.exponential_type(::Type{<:AbstractCliffordNumber})
    CliffordNumbers.exponential_type(x::AbstractCliffordNumber)

Returns the type expected when exponentiating a Clifford number. This is an `EvenCliffordNumber` if
the nonzero grades of the input are even, a `CliffordNumber` otherwise.
"""
@generated function exponential_type(::Type{C}) where {Q,C<:AbstractCliffordNumber{Q}}
    T = typeof(exp(zero(scalar_type(C))))
    if all(iseven, nonzero_grades(C))
        return :(EvenCliffordNumber{Q,$T,div(blade_count(Q), 2)})
    else
        return :(CliffordNumber{Q,$T,blade_count(Q)})
    end
end

exponential_type(x::AbstractCliffordNumber) = exponential_type(typeof(x))

# Odd grade Clifford numbers promote incorrectly due to broken type inference, see this issue:
# https://github.com/JuliaLang/julia/issues/53504
^(k::Union{KVector,OddCliffordNumber}, n::Integer) = convert(exponential_type(k), k)^n

"""
    CliffordNumbers.exp_taylor(x::AbstractCliffordNumber, order = Val(16))

Calculates the exponential of `x` using a Taylor expansion up to the specified order. In most cases,
12 is as sufficient number.

# Notes

16 iterations is currently used because the number of loop iterations is not currently a performance
bottleneck.
"""
@generated function exp_taylor(x::AbstractCliffordNumber, ::Val{N} = Val(16)) where N
    T = exponential_type(x)
    ex = quote
        # Initial term and result come from the zero exponent
        result = term = one($T)
        # Scale down the magnitude of x to prevent overflow
        # Use a power of 2 to simplify the later exponentiation
        # TODO: do this with abs2 instead of abs?
        r = div(exponent(2*abs2(x) + 1), 2)
        # Promote the argument to the final return type
        y = convert($T, x / (2^r))
    end
    # Unroll the iterations (16 seems to be enough for Float64, 12 for Float32)
    for n in 1:N
        ex = quote
            $ex
            # Use fast multiplication kernel directly
            y_scaled = y * convert(scalar_type($T), 1 // $n)
            term = mul(term, y_scaled)
            # Add the term from this loop to the final result
            result = result + term
        end
    end
    return quote
        $ex
        # Use the identity exp(x) = exp(x/s)^s
        for _ in 1:r
            result = mul(result, result)
        end
        return result
    end
end

# We can get away with fewer iterations for Float32 and Float16
exp_taylor(x::AbstractCliffordNumber{<:Any,Float32}) = exp_taylor(x, Val(12))
exp_taylor(x::AbstractCliffordNumber{<:Any,Float16}) = exp_taylor(x, Val(8))

"""
    exp(x::AbstractCliffordNumber{Q})

Returns the natural exponential of a Clifford number.

For special cases where m squares to a scalar, the following shortcuts can be used to calculate
`exp(x)`:
  * When x^2 < 0: `exp(x) === cos(abs(x)) + x * sin(abs(x)) / abs(x)`
  * When x^2 > 0: `exp(x) === cosh(abs(x)) + x * sinh(abs(x)) / abs(x)`
  * When x^2 == 0: `exp(x) == 1 + x`

See also: [`exppi`](@ref), [`exptau`](@ref).
"""
function exp(x::AbstractCliffordNumber)
    T = exponential_type(x)
    Tx = convert(T, x)
    sq = mul(Tx, Tx)
    if isscalar(sq)
        mag = abs(x)
        scalar(sq) < 0 && return T(cos(mag))  + Tx * (sin(mag) / mag)
        scalar(sq) > 0 && return T(cosh(mag)) + Tx * (sinh(mag) / mag)
        return one(T) + Tx
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
    Tx = convert(T, x)
    sq = mul(Tx, Tx)
    if isscalar(sq)
        mag = abs(x)
        scalar(sq) < 0 && return T(cospi(mag))  + Tx * (sinpi(mag) / mag)
        # TODO: is this really more accurate?
        scalar(sq) > 0 && return T(cosh(π*mag)) + Tx * (sinh(π*mag) / mag)
        return one(T) + Tx
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
