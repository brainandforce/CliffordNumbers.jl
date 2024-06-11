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

Base.inv(x::KVector{0,Q}) where Q = KVector{0,Q}(inv(scalar(x)))
Base.inv(x::KVector{1,Q}) where Q = versor_inverse(x)

/(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q = x * inv(y)
\(x::AbstractCliffordNumber{Q}, y::AbstractCliffordNumber{Q}) where Q = inv(x) * y
