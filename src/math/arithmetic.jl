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

#---Addition, negation, and subtraction------------------------------------------------------------#

+(x::T, y::T) where T<:AbstractCliffordNumber = T(map(+, Tuple(x), Tuple(y)))

function +(x::AbstractCliffordNumber, y::BaseNumber)
    T = promote_type(typeof(x), typeof(y))
    b = BitIndices(T)
    data = ntuple(i -> x[b[i]] + (iszero(grade(b[i])) * y), Val(nblades(T)))
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
    data = ntuple(i -> x[b[i]] - (iszero(grade(b[i])) * y), Val(nblades(T)))
    return T(data)
end

function -(x::BaseNumber, y::AbstractCliffordNumber)
    T = promote_type(typeof(x), typeof(y))
    b = BitIndices(T)
    data = ntuple(i -> (iszero(grade(b[i])) * x) - y[b[i]], Val(nblades(T)))
    return T(data)
end

# TODO: is it more efficient to define some more specific methods for some types?
