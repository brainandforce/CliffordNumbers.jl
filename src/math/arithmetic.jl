#---Equality and approximate equality--------------------------------------------------------------#

==(x::T, y::T) where T<:AbstractCliffordNumber = Tuple(x) == Tuple(y)
isequal(x::T, y::T) where T<:AbstractCliffordNumber = isequal(Tuple(x), Tuple(y))

# TODO: does this avoid unnecessary comparisons or construction?
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

#---Other properties-------------------------------------------------------------------------------#

# Already works: Base.isfinite(x::AbstractCliffordNumber) = all(isfinite, Tuple(x))
isinf(x::AbstractCliffordNumber) = any(isinf, Tuple(x))
isnan(x::AbstractCliffordNumber) = any(isnan, Tuple(x))

isreal(x::AbstractCliffordNumber) = isreal(scalar(x)) && isscalar(x)
isreal(x::AbstractCliffordNumber{<:Any,<:Real}) = isscalar(x)

isinteger(x::AbstractCliffordNumber) = isinteger(scalar(x)) && isscalar(x)
isinteger(x::AbstractCliffordNumber{<:Any,<:Integer}) = isscalar(x)

iseven(x::AbstractCliffordNumber) = isreal(x) && iseven(real(scalar(x)))
isodd(x::AbstractCliffordNumber) = isreal(x) && isodd(real(scalar(x)))

#---Addition, negation, and subtraction------------------------------------------------------------#

+(x::T, y::T) where T<:AbstractCliffordNumber = T(map(+, Tuple(x), Tuple(y)))

function -(x::AbstractCliffordNumber)
    data = (-).(Tuple(x))
    return similar_type(x, eltype(data))(data)
end

-(x::T, y::T) where T<:AbstractCliffordNumber = x + (-y)

# TODO: is it more efficient to define some more specific methods for some types?
