#---Efficient multiplication kernels---------------------------------------------------------------#
"""
    CliffordNumbers.bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}})
    CliffordNumbers.bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q})
    
    CliffordNumbers.bitindex_shuffle(B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q})
    CliffordNumbers.bitindex_shuffle(B::BitIndices{Q}, a::BitIndex{Q})

Performs the multiplication `-a * b` for each element of `B` for the above ordering, or `-b * a` for
the below ordering, generating a reordered `NTuple` of `BitIndex{Q}` objects suitable for
implementing a geometric product.
"""
@inline function bitindex_shuffle(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> a' * b, B)
end

@inline function bitindex_shuffle(B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q}
    return map(b -> b * a', B)
end

bitindex_shuffle(a::BitIndex{Q}, B::BitIndices{Q}) where Q = bitindex_shuffle(a, Tuple(B))
bitindex_shuffle(B::BitIndices{Q}, a::BitIndex{Q}) where Q = bitindex_shuffle(Tuple(B), a)

"""
    CliffordNumbers.nondegenerate_mask(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}})

Constructs a Boolean mask which is `false` for any multiplication that squares a degenerate blade;
`true` otherwise.
"""
function nondegenerate_mask(a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> nondegenerate_mult(a, b), B)
end

#=
"""
    CliffordNumbers.widen_grade_for_mul(x::AbstractCliffordNumber)

Widens `x` to an `EvenCliffordNumber`, `OddCliffordNumber`, or `CliffordNumber` as appropriate
for the fast multiplication kernel.
"""
widen_grade_for_mul(x::Union{CliffordNumber,Z2CliffordNumber}) = x
widen_grade_for_mul(k::KVector{K}) where K = Z2CliffordNumber{isodd(K)}(k)

# Generic fallback for future user-defined types
function widen_grade_for_mul(x::AbstractCliffordNumber)
    all(iseven, nonzero_grades(x)) && return EvenCliffordNumber(x)
    all(iodd, nonzero_grades(x)) && return OddCliffordNumber(x)
    return CliffordNumber(x)
end
=#

#---Grade filters----------------------------------------------------------------------------------#
"""
    CliffordNumbers.GradeFilter{S}

A type that can be used to filter certain products of blades in a geometric product multiplication.
The type parameter `S` must be a `Symbol`. The single instance of `GradeFilter{S}` is a callable
object which implements a function that takes two or more `BitIndex{Q}` objects `a` and `b` and
returns `false` if the product of the blades indexed is zero.

To implement a grade filter for a product function `f`, define the following method:
    (::GradeFilter{:f})(::BitIndex{Q}, ::BitIndex{Q})
    # Or if the definition allows for more arguments
    (::GradeFilter{:f})(::BitIndex{Q}...) where Q
"""
struct GradeFilter{S}
    GradeFilter{S}() where S = (@assert S isa Symbol "Type parameter must be a Symbol."; new())
end

(::GradeFilter{S})(args...) where S = error("This filter has not been implemented.")

(::GradeFilter{:*})(args::BitIndex{Q}...) where Q = true

(::GradeFilter{:∧})(args::BitIndex{Q}...) where Q = has_wedge(args...)

(::GradeFilter{:⨼})(a::BitIndex{Q}, b::BitIndex{Q}) where Q = (grade(b) - grade(a)) == grade(a*b)
(::GradeFilter{:⨽})(a::BitIndex{Q}, b::BitIndex{Q}) where Q = (grade(a) - grade(b)) == grade(a*b)

function (::GradeFilter{:dot})(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    return abs(grade(a) - grade(b)) == grade(a*b)
end

const ContractionGradeFilters = Union{GradeFilter{:⨼},GradeFilter{:⨽},GradeFilter{:dot}}

"""
    CliffordNumbers.mul_mask(F::GradeFilter, a::BitIndex{Q}, B::NTuple{L,BitIndices{Q}})
    CliffordNumbers.mul_mask(F::GradeFilter, B::NTuple{L,BitIndices{Q}}, a::BitIndex{Q})

    CliffordNumbers.mul_mask(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q})
    CliffordNumbers.mul_mask(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q})

Generates a `NTuple{L,Bool}` which is `true` whenever the multiplication of the blade indexed by `a`
and blades indexed by `B` is nonzero. `false` is returned if the grades multiply to zero due to the
squaring of a degenerate component, or if they are filtered by `F`.
"""
function mul_mask(F::GradeFilter, a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return map(b -> F(a,b) & nondegenerate_mult(a,b), B)
end

function mul_mask(F::GradeFilter, B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q}
    return map(b -> F(b,a) & nondegenerate_mult(b,a), B)
end

mul_mask(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q}) where Q = mul_mask(F, a, Tuple(B))
mul_mask(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q}) where Q = mul_mask(F, Tuple(B), a)

"""
    CliffordNumbers.mul_signs(F::GradeFilter, a::BitIndex{Q}, B::NTuple{L,BitIndices{Q}})
    CliffordNumbers.mul_signs(F::GradeFilter, B::NTuple{L,BitIndices{Q}}, a::BitIndex{Q})

    CliffordNumbers.mul_signs(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q})
    CliffordNumbers.mul_signs(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q})

Generates an `NTuple{L,Int8}` which represents the sign associated with the multiplication needed to
calculate components of a multiplication result.

This is equivalent to `sign.(B)` unless `F === CliffordNumbers.GradeFilter{:dot}()`.
"""
mul_signs(::GradeFilter, ::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q} = sign.(B)
mul_signs(::GradeFilter, B::NTuple{L,BitIndex{Q}}, ::BitIndex{Q}) where {L,Q} = sign.(B)

function mul_signs(::GradeFilter{:dot}, a::BitIndex{Q}, B::NTuple{L,BitIndex{Q}}) where {L,Q}
    return sign.(B) .* Int8(-1).^(grade.(B) .* (grade(a) .- grade.(B)))
end

function mul_signs(::GradeFilter{:dot}, B::NTuple{L,BitIndex{Q}}, a::BitIndex{Q}) where {L,Q}
    return sign.(B) .* Int8(-1).^(grade(a) .* (grade.(B) .- grade(a)))
end

mul_signs(F::GradeFilter, a::BitIndex{Q}, B::BitIndices{Q}) where Q = mul_signs(F, a, Tuple(B))
mul_signs(F::GradeFilter, B::BitIndices{Q}, a::BitIndex{Q}) where Q = mul_signs(F, Tuple(B), a)

#---Product return types---------------------------------------------------------------------------#
"""
    CliffordNumbers.product_return_type(::Type{X}, ::Type{Y}, [::GradeFilter{S}])

Returns a suitable type for representing the product of Clifford numbers of types `X` and `Y`. The
`GradeFilter{S}` argument allows for the return type to be changed depending on the type of product.
Without specialization on `S`, a type suitable for the geometric product is returned.
"""
@generated function product_return_type(
    ::Type{C1},
    ::Type{C2},
    ::GradeFilter{<:Any}
) where {Q,C1<:AbstractCliffordNumber{Q},C2<:AbstractCliffordNumber{Q}}
    c1_odd = all(isodd, nonzero_grades(C1))
    c2_odd = all(isodd, nonzero_grades(C2))
    c1_even = all(iseven, nonzero_grades(C1))
    c2_even = all(iseven, nonzero_grades(C2))
    # Parity: true for odd multivectors, false for even multivectors
    P = (c1_odd && c2_even) || (c1_even && c2_odd)
    T = promote_scalar_type(C1,C2)
    if (!c1_odd && !c1_even) || (!c2_odd && !c2_even)
        return :(CliffordNumber{Q})
    else
        return :(Z2CliffordNumber{$P,Q})
    end
end

# Extra implementation for k-vectors: special handling scalar and pseudoscalar arguments
# TODO: can we integrate this into the above function?
@generated function product_return_type(
    ::Type{C1},
    ::Type{C2},
    ::GradeFilter{<:Any}
) where {K1,K2,Q,C1<:KVector{K1,Q},C2<:KVector{K2,Q}}
    D = dimension(Q)
    # Handle the scalar and pseudoscalar cases
    if isone(nblades(C1))
        iszero(K1) && return :(KVector{K2,Q})
        K1 == D && return :(KVector{$D-K2,Q})
    elseif isone(nblades(C2))
        iszero(K2) && return :(KVector{K1,Q})
        K2 == D && return :(KVector{$D-K1,Q})
    end
    # Fall back to a Z2CliffordNumber with the right parity
    P = isodd(xor(K1, K2))
    return :(Z2CliffordNumber{$P,Q})
end

function product_return_type(
    ::Type{<:KVector{K1,Q}},
    ::Type{<:KVector{K2,Q}},
    ::GradeFilter{:∧}
) where {Q,K1,K2}
    # TODO: do we need this minimizing behavior?
    # Maybe we can let K exceed the expected value and return a zero-element multivector.
    K = min(K1 + K2, dimension(Q))
    return KVector{K,Q}
end

function product_return_type(
    ::Type{<:KVector{K1,Q}},
    ::Type{<:KVector{K2,Q}},
    ::ContractionGradeFilters
) where {Q,K1,K2}
    K = abs(K1 - K2)
    return KVector{K,Q}
end

function product_return_type(
    x::AbstractCliffordNumber,
    y::AbstractCliffordNumber,
    F::GradeFilter = GradeFilter{:*}()
)
    return product_return_type(typeof(x), typeof(y), F)
end

#---Geometric product------------------------------------------------------------------------------#
"""
    CliffordNumbers.mul(
        x::AbstractCliffordNumber{Q},
        y::AbstractCliffordNumber{Q},
        [F::GradeFilter = GradeFilter{:*}()]
    )

A fast geometric product implementation using generated functions for specific cases, and generic
methods which either convert the arguments or fall back to other methods.

The arguments to this function should all agree in scalar type `T`. The `*` function, which exposes
the fast geometric product implementation, promotes the scalar types of the arguments before
utilizing this kernel. The scalar multiplication operations are implemented using [`muladd`](@ref), 
allowing for hardware fma operations to be used when available.

The `GradeFilter` `F` allows for some blade multiplications to be excluded if they meet certain
criteria. This is useful for implementing products besides the geometric product, such as the wedge
product, which excludes multiplications between blades with shared vectors. Without a filter, this
kernel just returns the geometric product.
"""
@generated function mul(
    x::AbstractCliffordNumber{Q,T},
    y::AbstractCliffordNumber{Q,T},
    F::GradeFilter = GradeFilter{:*}()
) where {Q,T<:BaseNumber}
    C = product_return_type(x, y, F())
    ex = :($(zero_tuple(C)))
    # TODO: specialize kernel to favor iterating through the indices of the smaller argument
    for a in BitIndices(x)
        # Permute the indices of y so that the result coefficients are correctly ordered
        inds = bitindex_shuffle(a, BitIndices(C))
        # Filter out multiplications which necessarily go to zero
        x_mask = mul_mask(F(), a, inds)
        # Filter out indexing operations that automatically go to zero
        # This must be done manually since we want to work directly with tuples
        y_mask = map(in, grade.(inds), ntuple(Returns(nonzero_grades(y)), Val(nblades(C))))
        mask = x_mask .& y_mask
        # Don't append operations that won't actually do anything
        if any(mask)
            # Resolve BitIndex to an integer here to avoid having to call Base.to_index at runtime
            # This function cannot be inlined or unrolled for KVector arguments
            # But all values are known at compile time, so interpolate them into expressions
            ia = to_index(x, a)
            tuple_inds = to_index.(y, inds)
            signs = mul_signs(F(), a, inds)
            # Construct the tuples that contribute to the product
            x_tuple_ex = :(Tuple(x)[$ia] .* $(signs .* mask))
            y_tuple_ex = :(getindex.(tuple(Tuple(y)), $tuple_inds))
            # Combine the tuples using muladd operations
            ex = :(map(muladd, $x_tuple_ex, $y_tuple_ex, $ex))
        end
    end
    return :(($C)($ex))
end

#= Update (2024-05-07)
    The performance bottlenecks for KVector arguments have been (mostly) resolved.
    This was done by resolving all BitIndex objects to tuple indices at compile time, avoiding
    all calls to Base.to_index(::KVector, ::Int) at runtime (since this function cannot be unrolled
    or inlined at compile time).

    There might be some lingering performance issues with the order of differently sized arguments.
=#

function mul(
    x::AbstractCliffordNumber{Q},
    y::AbstractCliffordNumber{Q},
    F::GradeFilter = GradeFilter{:*}()
) where Q
    return @inline mul(scalar_promote(x, y)..., F)
end

function mul(
    x::AbstractCliffordNumber,
    y::AbstractCliffordNumber,
    ::GradeFilter{S} = GradeFilter{:mul}()
) where S
    throw(AlgebraMismatch(S, (x, y)))
end
