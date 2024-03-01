# Default conversion should check for exact representability
function convert(::Type{T}, x::AbstractCliffordNumber{Q}) where {Q,T<:AbstractCliffordNumber{Q}}
    has_grades_of(x, T) && return T(x)::T
    throw(InexactError(:convert, T, x))
end

function convert(::Type{T}, x::AbstractCliffordNumber) where T<:BaseNumber
    return isscalar(x) ? T(scalar(x)) : throw(InexactError(:convert, T, x))
end

#---Specialized conversion methods for certain representations and signatures----------------------#

function convert(::Type{T}, x::AbstractCliffordNumber{QFComplex,<:Real}) where T<:BaseNumber
    return convert(T, x[scalar_index(x)] + x[pseudoscalar_index(x)] * im)
end

function convert(::Type{T}, z::Complex) where T<:AbstractCliffordNumber{QFComplex,<:Real}
    return convert(T, CliffordNumber{QFComplex}(real(z), imag(z)))
end

# k-vectors of grade 0 are scalars
convert(::Type{T}, k::KVector{0}) where T<:BaseNumber = convert(T, only(k.data))
