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
    # TODO: do we need to account for the sizes of the inputs?
    if isfinite(d)
        # The norm of the difference should not exceed the tolerance
        return d <= max(atol, rtol*max(norm(x), norm(y)))
    else
        # Compare each element of each Clifford number
        (tx, ty) = Tuple.(promote(x, y))
        return all(isapprox.(tx, ty; atol, rtol, nans, norm))
    end
    #= TODO:
    If we strip away the `isfinite` if block to simplify the function, we get *worse* performance!
    This is probably due to the closure bug, but why does the if block avoid the problem?
    =#
end

isapprox(x::AbstractCliffordNumber, y::Number; kwargs...) = isapprox(promote(x, y)...; kwargs...)
isapprox(x::Number, y::AbstractCliffordNumber; kwargs...) = isapprox(promote(y, x)...; kwargs...)

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

+(x::BaseNumber, y::AbstractCliffordNumber) = y + x

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
            data = Tuple(y') .* $op(x, sum(Tuple(y) .* Tuple(y')))
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

*(x::AbstractCliffordNumber, y::AbstractCliffordNumber) = throw(AlgebraMismatch(*, (x, y)))

# KVector{0,Q} is just a scalar compatible with an AbstractCliffordNumber{Q}
# TODO: this may be best in an eval block with all other products
*(k::KVector{0,Q}, x::AbstractCliffordNumber{Q}) where Q = only(Tuple(k)) * x
*(x::AbstractCliffordNumber{Q}, k::KVector{0,Q}) where Q = x * only(Tuple(k))
*(k::KVector{0,Q}, l::KVector{0,Q}) where Q = KVector{0,Q}(only(Tuple(k)) * only(Tuple(l)))

@inline (x::AbstractCliffordNumber)(y::AbstractCliffordNumber) = x * y

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
    abs2(x::AbstractCliffordNumber{Q,T}) -> Union{Real,Complex}

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
    CliffordNumbers.dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the dot product of `x` and `y`.

For basis blades `A` of grade `m` and `B` of grade `n`, the dot product is equal to the left
contraction when `m >= n` and is equal to the right contraction (up to sign) when `n >= m`.

# Why is this function not exported?

The LinearAlgebra package also defines a `dot` function, and if both packages are used together,
this will cause a name conflict if `CliffordNumbers.dot` is exported. In the future, we will try to
resolve this without requiring a LinearAlgebra dependency.

Additionally, there is reason to prefer the use of the left and right contractions over the dot
product because the contractions require fewer exceptions in their definitions and properties.
"""
function dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q 
    return mul(scalar_promote(x, y)..., GradeFilter{:dot}())
end

const left_contraction = ⨼
const right_contraction = ⨽

"""
    CliffordNumbers.hestenes_dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Returns the Hestenes product: this is equal to the dot product given by `dot(x, y)` but is equal to
to zero when either `x` or `y` is a scalar.

# Why is this function not exported?

In almost every case, left and right contractions are preferable - the dot product and the Hestenes
product are less regular in algebraic sense. It is provided for the sake of exact reproducibility of
results which use it.
"""
function hestenes_dot(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
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
    ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    commutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the commutator (or antisymmetric) product, equal to `1//2 * (x*y - y*x)`.

Note that the commutator product, in general, is *not* equal to the wedge product, which may be
invoked with the `wedge` function or the `∧` operator.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function ×(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y - y*x)
end

const commutator = ×

"""
    ⨰(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})
    anticommutator(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q})

Calculates the anticommutator (or symmetric) product, equal to `1//2 * (x*y + y*x)`.

Note that the dot product, in general, is *not* equal to the anticommutator product, which may be
invoked with `dot`. In some cases, the preferred operators might be the left and right contractions,
which use infix operators `⨼` and `⨽` respectively.

# Type promotion

Because of the rational `1//2` factor in the product, inputs with scalar types subtyping `Integer`
will be promoted to `Rational` subtypes.
"""
function ⨰(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q
    return 1//2 * (x*y + y*x)
end

const anticommutator = ⨰

#---Inverses and division--------------------------------------------------------------------------#
"""
    CliffordNumbers.InverseException(msg)

No inverse exists for the given input. Optional parameter `msg` is a descriptive error string.
"""
struct InverseException <: Exception
    msg::String
end

function Base.showerror(io::IO, ex::InverseException)
    print(io, typeof(ex), ": ", ex.msg)
end

"""
    CliffordNumbers.versor_inverse(x::AbstractCliffordNumber)

Calculates the versor inverse of `x`, equal to `x' / abs2(x)`.

The versor inverse is only guaranteed to be an inverse for versors. Not all Clifford
numbers have a well-defined inverse, (for instance, in algebras with 2 or more positive-squaring,
dimensions, 1 + e₁ has no inverse). To validate the result, use `inv(x)` instead.
"""
versor_inverse(x::AbstractCliffordNumber) = x' / abs2(x)

"""
    inv(x::AbstractCliffordNumber) -> AbstractCliffordNumber

Calculates the inverse of `x`, if it exists, using the versor inverse formula `x' / abs2(x)`. The
result is tested to check if its left and right products with `x` are approximately 1, and a 
`CliffordNumbers.InverseException` is thrown if this test does not pass.
"""
function Base.inv(x::AbstractCliffordNumber)
    invx = versor_inverse(x)
    # Pseudovectors and pseudoscalars always have inverses
    # Scalars and vectors are handled with separate methods below
    if x isa KVector && dimension(signature(x)) - grade(x) in (0,1)
        return invx
    # Explicitly check that the inverse exists
    elseif x * invx ≈ 1 && invx * x ≈ 1
        return invx
    else
        throw(
            InverseException(
                "The result of x' / abs2(x) was not an inverse.\n" * 
                "Note that not all Clifford numbers will have an inverse."
            )
        )
    end
end

Base.inv(x::KVector{0,Q}) where Q = KVector{0,Q}(inv(only(Tuple(x))))
Base.inv(x::KVector{1,Q}) where Q = versor_inverse(x)

/(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q = x * inv(y)
\(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q = inv(x) * y

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
