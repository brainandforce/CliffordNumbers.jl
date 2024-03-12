#---Efficient multiplication kernels---------------------------------------------------------------#

@inline bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q} = map(b -> -a*b, B)
bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q}) where Q = map(b -> a*b, Tuple(B))

function _ndmult(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> nondegenerate_mult(a, b), B)
end

@generated function mul(
    ::Type{C},
    x::T,
    y::T,
) where {Q,C<:AbstractCliffordNumber{Q},T<:AbstractCliffordNumber{Q}}
    BC = Tuple(BitIndices(C))
    z = zero_tuple(C)
    ex = :($z)
    for a in BitIndices(x)
        inds = bitindex_shuffle(a, BC)
        nd = _ndmult(a, BC)
        ex = :(map(muladd, x[$a] .* $nd, y[$inds], $ex))
    end
    return :(C($ex))
end
