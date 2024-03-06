# Default conversion should check for exact representability
function convert(T::Type{<:AbstractCliffordNumber}, x::AbstractCliffordNumber)
    result = T(x)::T
    return (has_grades_of(x, result) ? result : throw(InexactError(:convert, T, x)))
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
