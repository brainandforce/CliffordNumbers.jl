"""
    KVector{K,Q,T,L} <: AbstractCliffordNumber{Q,T}

A multivector consisting only linear combinations of basis blades of grade `K` - in other words,
a k-vector.

k-vectors have `binomial(dimension(Q), K)` components.
"""
struct KVector{K,Q,T<:BaseNumber,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function KVector{K,Q,T,L}(x::Tuple) where {K,Q,T,L}
        @assert 0 <= K <= dimension(Q) "K can only range from 0 to $(dimension(Q)) (got $K)."
        check_element_count(binomial(dimension(Q), K), L, x)
        return new{K,Q,T,L}(x)
    end
end

#---Constructors-----------------------------------------------------------------------------------#

KVector{K,Q,T}(x::Tuple{Vararg{BaseNumber,L}}) where {K,Q,T,L} = KVector{K,Q,T,L}(x)
KVector{K,Q}(x::Tuple{Vararg{T}}) where {K,Q,T<:BaseNumber} = KVector{K,Q,T}(x)

# Automatically convert arguments to a common type
KVector{K,Q}(x::Tuple{Vararg{BaseNumber}}) where {K,Q} = KVector{K,Q}(promote(x...))

#---Number of elements-----------------------------------------------------------------------------#

length(::Type{KVector{K,Q,T,L}}) where {K,Q,T,L} = L
length(::Type{<:KVector{K,Q}}) where {K,Q} = binomial(dimension(Q), K)
length(x::KVector) = length(typeof(x))

#---Indexing---------------------------------------------------------------------------------------#

size(::BitIndices{Q,<:KVector{K}}) where {Q,K} = tuple(length(KVector{K,Q}))

"""
    grade(::Type{<:KVector{K}}) = K
    grade(x::KVector{K}) = k

Returns the grade represented by a `KVector{K}`, which is K.
"""
grade(::Type{<:KVector{K}}) where K = K
grade(x::KVector) = grade(typeof(x))

nonzero_grades(::Type{<:KVector{K}}) where K = K:K

@inline function getindex(b::BitIndices{Q,<:KVector{K}}, i::Integer) where {Q,K}
    @boundscheck checkbounds(b, i)
    return BitIndex{Q}(signbit(i-1), unsigned(hamming_number(K, i)))
end

@inline function to_index(C::Type{<:KVector{K,Q}}, b::BitIndex{Q}) where {K,Q}
    # Default to 1 as a valid index for any KVector instance
    i = 1
    for n in 1:length(C)
        is_same_blade(b, (@inbounds BitIndices(C)[n])) && (i = n)
    end
    return i
end

@inline to_index(k::KVector{K,Q}, b::BitIndex{Q}) where {K,Q} = to_index(typeof(k), b)

@inline function getindex(k::KVector{K,Q}, b::BitIndex{Q}) where {K,Q}
    return (@inbounds k.data[to_index(k, b)]) * sign(b) * (grade(b) === K)
end

#---Multiplicative identity and pseudoscalar-------------------------------------------------------#

one(K::Type{<:KVector{<:Any,Q}}) where Q = KVector{0,Q}(scalar_type(K)(true))
# Default implementation for AbstractCliffordNumber subtypes if not explicitly given
one(C::Type{<:AbstractCliffordNumber{Q}}) where Q = one(KVector{0,Q,scalar_type(C)})

"""
    pseudoscalar(C::Type{<:AbstractCliffordNumber{Q,T}}) -> KVector{dimension{Q},Q,T}
    pseudoscalar(x::AbstractCliffordNumber{Q})

Returns the pseudoscalar associated with the signature `Q` of the argument. The result will have the
same scalar type as the input if such information is given; otherwise it will use `Bool` as the
scalar type so as to promote to any other numeric type.
"""
function pseudoscalar(C::Type{<:AbstractCliffordNumber{Q}}) where Q
    return KVector{dimension(Q),Q}(scalar_type(C)(true))
end

pseudoscalar(x::AbstractCliffordNumber) = pseudoscalar(typeof(x))

#---Similar types----------------------------------------------------------------------------------#

function similar_type(::Type{<:KVector{K}}, ::Type{T}, ::Val{Q}) where {K,Q,T<:BaseNumber}
    return KVector{K,Q,T,binomial(dimension(Q),K)}
end

function similar_type(::Type{<:KVector}, ::Type{T}, ::Val{Q}, ::Val{K}) where {K,Q,T<:BaseNumber}
    return KVector{K,Q,T,binomial(dimension(Q),K)}
end

#---Show methods-----------------------------------------------------------------------------------#

short_typename(::Type{<:KVector{K,Q,T}}) where {K,Q,T} = KVector{K,Q,T}

#---Special Z2CliffordNumber constructor-----------------------------------------------------------#

# Automatically infer if we want an even or odd Clifford number
Z2CliffordNumber(x::KVector{K}) where K = Z2CliffordNumber{isodd(K)}(x)
