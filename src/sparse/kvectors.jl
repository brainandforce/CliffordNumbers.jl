"""
    KVector{K,Q,T,L} <: AbstractCliffordNumber{Q,T}

A multivector consisting only linear combinations of basis blades of grade `K` - in other words,
a k-vector.

k-vectors have `binomial(dimension(Q), K)` components.
"""
struct KVector{K,Q,T,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function KVector{K,Q,T,L}(x) where {K,Q,T,L}
        sz = binomial(dimension(Q), K)
        @assert length(x) == L == sz string(
            "Incorrect number of components: ", K, "-vectors of ", Q, " have ", sz, " components."
        )
        return new{K,Q,T,L}(x)
    end
end

KVector{K,Q,T}(x::NTuple{L,<:BaseNumber}) where {K,Q,T,L} = KVector{K,Q,T,L}(x)
KVector{K,Q,T}(x::Vararg{BaseNumber,L}) where {K,Q,T,L} = KVector{K,Q,T,L}(x)

KVector{K,Q}(x::NTuple{L,T}) where {K,Q,T,L} = KVector{K,Q,T,L}(x)
KVector{K,Q}(x::Vararg{T,L}) where {K,Q,T,L} = KVector{K,Q,T,L}(x)

#---Number of elements-----------------------------------------------------------------------------#
import Base: length

length(::Type{KVector{K,Q,T,L}}) where {K,Q,T,L} = L
length(::Type{<:KVector{K,Q}}) where {K,Q} = binomial(dimension(Q), K)
length(x::KVector) = length(typeof(x))

#---Indexing---------------------------------------------------------------------------------------#

Base.size(::BitIndices{Q,<:KVector{K}}) where {Q,K} = tuple(length(KVector{K,Q}))

function Base.getindex(b::BitIndices{Q,<:KVector{K}}, i::Integer) where {Q,K}
    @boundscheck checkbounds(b, i)
    return BitIndex{Q}(signbit(i-1), unsigned(hamming_number(K, i)))
end

function Base.getindex(k::KVector{K,Q}, b::BitIndex{Q}) where {K,Q}
    # Indices with mismatched grades are always zero
    grade(b) === K || return zero(eltype(k))
    return k.data[findfirst(isequal(abs(b)), BitIndices(typeof(k)))] * sign(b)
end

#---Zero elements----------------------------------------------------------------------------------#
import Base: zero

zero(::Type{KVector{K,Q,T,L}}) where {K,Q,T,L} = KVector{K,Q,T,L}(NTuple{L,T}(zero(T) for n in 1:L))
zero(S::Type{KVector{K,Q,T}}) where {K,Q,T} = zero(S{length(S)})
zero(::Type{KVector{K,Q}}) where {K,Q} = zero(KVector{K,Q,Bool,length(KVector{K,Q})})
zero(x::KVector) = zero(typeof(x))

#---Show methods-----------------------------------------------------------------------------------#
import Base: show, summary

function show(io::IO, k::KVector{K}) where K
    print(io, "KVector{", K, ",", algebra(k), ",", eltype(k), "}", k.data)
end

function summary(io::IO, k::KVector{K}) where K
    println(io, "KVector{", K, ",", algebra(k), ",", eltype(k), "}:")
end
