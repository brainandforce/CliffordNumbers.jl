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

# If the exponent is a constant, we can leverage that information to promote conservatively
@inline function Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{n}) where n
    # Handle negative cases: invert x, flip the sign of n, and retry
    n < 0 && return Base.literal_pow(^, inv(x), Val(abs(n)))
    # Handle positive cases: force promotion to even multivectors
    n > 0 && return (x*x)^(div(n,2)) * Base.literal_pow(^, x, Val(rem(n,2)))
    return one(x)
end

# Overload Base.literal_pow for common cases
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{2}) = x*x

@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{-1}) = inv(x)
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{-2}) = (i = inv(x); i*i)

# It appears that exponentiation with `Bool` does not get converted to Base.literal_pow
# But they're defined here anyway
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{false}) = one(x)
@inline Base.literal_pow(::typeof(^), x::AbstractCliffordNumber, ::Val{true}) = x

# Odd grade Clifford numbers promote incorrectly by default, because typeof(one(x)) != typeof(x)
# See this issue: https://github.com/JuliaLang/julia/issues/53504
^(x::C, n::Integer) where C<:Union{KVector,OddCliffordNumber} = convert(exponential_type(C), x)^n

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
