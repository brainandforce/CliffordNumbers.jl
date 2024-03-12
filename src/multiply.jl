#---Efficient multiplication kernels---------------------------------------------------------------#

@inline bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q} = map(b -> -a*b, B)
bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q}) where Q = map(b -> a*b, Tuple(B))

function _ndmult(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> nondegenerate_mult(a, b), B)
end

#=
    Specialized kernel for multiplying two CliffordNumber or Z2CliffordNumber instances
    
    NOTE: this is much faster if the arguments are the same type
    (CliffordNumber or Z2CliffordNumber; parity is irrelevant for the latter)
    TODO: handle the case where the smaller argument comes first

    NOTE: this is much faster if the first input is smaller than the second input
=#
@generated function mul(
    ::Type{C},
    x::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}},
    y::Union{CliffordNumber{Q,T},Z2CliffordNumber{<:Any,Q,T}}
) where {Q,T,C<:AbstractCliffordNumber{Q,T}}
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
