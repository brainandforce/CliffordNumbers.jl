"""
    CliffordNumbers.signmask([T::Type{<:Integer} = UInt], [signbit::Bool = true]) -> T

Generates a signmask, or a string of bits where the only 1 bit is the sign bit. If `signbit` is set
to false, this returns zero (or whatever value is represented by all bits being 0).
"""
@inline signmask(T::Type{<:Integer}, signbit::Bool = true) = bitreverse(T(signbit))
@inline signmask(x::Integer, signbit::Bool = true) = bitreverse(typeof(x)(signbit))

"""
    BitIndex{Q}

A representation of an index corresponding to a basis blade of the geometric algebra with quadratic
form `Q`.
"""
struct BitIndex{Q}
    i::UInt
    BitIndex{Q}(i::Unsigned) where Q = new(UInt(i % blade_count(Q)) | (UInt(i) & signmask(i)))
end

# Construct with the sign bit as a separate argument
@inline function BitIndex{Q}(signbit::Bool, blade::Unsigned) where Q
    return BitIndex{Q}(UInt(blade % blade_count(Q)) | signmask(UInt, signbit))
end

# UInt(i::BitIndex) strips the sign info. Use i.i to directly access the internal field
Base.UInt(i::BitIndex) = i.i & ~signmask(Int)
Base.Int(i::BitIndex) = Int(UInt(i))

#---Treat as scalar for the sake of broadcasting---------------------------------------------------#

Broadcast.broadcastable(b::BitIndex) = tuple(b)

#---Convenience constructors-----------------------------------------------------------------------#
"""
    CliffordNumbers._sort_with_parity!(v::AbstractVector{<:Real}) -> Tuple{typeof(v),Bool}

Performs a parity-tracking insertion sort of `v`, which modifies `v` in place. The function returns
a tuple containing `v` and the parity, which is `true` for an odd permutation and `false` for an
even permutation. This is implemented with a modified insertion sort algorithm.
"""
function _sort_with_parity!(v::AbstractVector{<:Real})
    swaps = 0
    # Algorithm adapted from Julia Base.Sort._sort!
    for i in axes(v,1)[2:end]
        j = i
        x = v[i]
        # j will decrement until it equals the first index of v
        while j > firstindex(v)
            # Look behind element i
            y = v[j-1]
            # It's sorted when x is greater than or equal to y
            x < y || break
            # Insert y where x was
            v[j] = y
            # Move to the previous index for comparison
            j -= 1
            # Track number of swaps
            swaps += 1
        end
        # Complete the swap of x and y values
        v[j] = x
    end
    return (v, isodd(swaps))
end

function _sort_with_parity(t::NTuple{L}) where L
    (v, sign) = _sort_with_parity!(collect(t))
    return (ntuple(i -> v[i], Val(L)), sign)
end

function _bitindex(::Val{S}, t::NTuple) where S
    @assert all(in(eachindex(S)), t) "Some of the indices are not valid for the given signature."
    (t, parity) = _sort_with_parity(t)
    i = signmask(UInt, parity)
    for x in t
        i = xor(i, <<(1, x - firstindex(S)))
    end
    return BitIndex{S}(i)
end

"""
    BitIndex(::Val{Q}, i::Integer...)
    BitIndex(x, i::Integer...) = BitIndex(signature(x), i...)

Constructs a `BitIndex{Q}` from a list of integers that represent the basis 1-vectors of the space.
For objects `x` that are not metric signatures, `Q` is equal to `signature(x)`.

This package uses a lexicographic convention for basis blades: in the algebra of physical space, the
basis bivectors are {e₁e₂, e₁e₃, e₂e₃}. The sign of the `BitIndex{Q}` is negative when the parity of
the basis vector permutation is odd.
"""
BitIndex(S::Val, i::Integer...) = _bitindex(S, promote(i...))
BitIndex(x, i::Integer...) = BitIndex(Val(signature(x)), i...)

#---Show method------------------------------------------------------------------------------------#

function show(io::IO, b::BitIndex{Q}) where Q
    print(io, "-"^signbit(b), BitIndex, "(Val(", Q, ")", iszero(Int(b)) ? ")" : ", ")
    iszero(UInt(b)) && return nothing
    found_first_vector = false
    for a in 1:min(dimension(Q), 8*sizeof(UInt) - 1)
        if !iszero(UInt(b) & <<(1, a-1))
            found_first_vector && print(io, ", ")
            print(io, eachindex(Q)[a])
            found_first_vector = true
        end
    end
    print(io, ")")
end

#---Other useful functions-------------------------------------------------------------------------#

signature(::BitIndex{Q}) where Q = Q

signbit(i::BitIndex) = !iszero(i.i & signmask(UInt))
sign(i::BitIndex) = Int8(-1)^signbit(i)

-(i::BitIndex) = typeof(i)(!signbit(i), UInt(i))
abs(i::BitIndex) = typeof(i)(UInt(i))

"""
    grade(i::BitIndex) -> Int

Returns the grade of the basis blade represented by `i`, which ranges from 0 to the dimension of the
space represented by `i` (equal to `dimension(signature(i))`).
"""
grade(i::BitIndex) = count_ones(UInt(i))

"""
    CliffordNumbers.is_same_blade(a::BitIndex{Q}, b::BitIndex{Q})

Checks if `a` and `b` perform identical indexing up to sign.
"""
is_same_blade(a::T, b::T) where T<:BitIndex = (a.i << 1) === (b.i << 1)
is_same_blade(a::BitIndex, b::BitIndex) = false

"""
    scalar_index(x::AbstractCliffordNumber{Q}) -> BitIndex{Q}

Constructs the `BitIndex` used to obtain the scalar (grade zero) portion of `x`.
"""
scalar_index(::Type{<:AbstractCliffordNumber{Q}}) where Q = BitIndex{Q}(false, UInt(0))
scalar_index(x::AbstractCliffordNumber) = scalar_index(typeof(x))

"""
    pseudoscalar_index(x::AbstractCliffordNumber{Q}) -> BitIndex{Q}

Constructs the `BitIndex` used to obtain the pseudoscalar (highest grade) portion of `x`.
"""
pseudoscalar_index(::Type{<:AbstractCliffordNumber{Q}}) where Q = BitIndex{Q}(false, typemax(UInt))
pseudoscalar_index(x::AbstractCliffordNumber) = pseudoscalar_index(typeof(x))

#---Grade dependent sign inversion-----------------------------------------------------------------#
"""
    adjoint(i::BitIndex) = reverse(i::BitIndex) = i' -> BitIndex
    adjoint(x::AbstractCliffordNumber) = reverse(x::AbstractCliffordNumber) = x' -> typeof(x)

Performs the reverse operation on the basis blade indexed by `b` or the Clifford number `x`. The 
sign of the reverse depends on the grade of the basis blade `g`, and is positive for `g % 4 in 0:1`
and negative for `g % 4 in 2:3`.
"""
adjoint(i::BitIndex) = typeof(i)(xor(signbit(i), !iszero(grade(i) & 2)), UInt(i))
reverse(i::BitIndex) = adjoint(i)
@doc (@doc adjoint(::BitIndex)) reverse(::BitIndex)

"""
    grade_involution(i::BitIndex) -> BitIndex
    grade_involution(x::AbstractCliffordNumber) -> typeof(x)

Calculates the grade involution of the basis blade indexed by `b` or the Clifford number `x`. This
effectively reflects all of the basis vectors of the space along their own mirror operation, which
makes elements of odd grade flip sign.
"""
grade_involution(i::BitIndex) = typeof(i)(xor(signbit(i), isodd(grade(i))), UInt(i))

"""
    conj(i::BitIndex) -> BitIndex
    conj(x::AbstractCliffordNumber) -> typeof(x)

Calculates the Clifford conjugate of the basis blade indexed by `b` or the Clifford number `x`. This
is equal to `grade_involution(reverse(x))`.
"""
conj(i::BitIndex) = typeof(i)(xor(signbit(i), !iszero((grade(i) + 1) & 2)), UInt(i))

#---Multiplication tools---------------------------------------------------------------------------#

# NOTE: <<(UInt(1), x) is faster/easier to inline than UInt(2)^x

function positive_square_bits(S::Metrics.AbstractSignature)
    return sum((S[x] === Int8(+1)) * <<(UInt(1), x - firstindex(S)) for x in eachindex(S))
end

function negative_square_bits(S::Metrics.AbstractSignature)
    return sum((S[x] === Int8(-1)) * <<(UInt(1), x - firstindex(S)) for x in eachindex(S))
end

function zero_square_bits(S::Metrics.AbstractSignature)
    return sum((S[x] === Int8(0)) * <<(UInt(1), x - firstindex(S)) for x in eachindex(S))
end

"""
    CliffordNumbers.signbit_of_square(b::BitIndex) -> Bool

Returns the signbit associated with squaring the basis blade indexed by `b`.

For basis blades squaring to zero, the result is not meaningful.
"""
@inline function signbit_of_square(b::BitIndex{Q}) where Q
    return xor(!iszero(grade(b) & 2), isodd(count_ones(UInt(b) & negative_square_bits(Q))))
end

"""
    CliffordNumbers.nondegenerate_square(b::BitIndex) -> Bool

Returns `false` if squaring the basis blade `b` is zero due to a degenerate component, `true` 
otherwise. For a nondegenerate metric, this is always `true`.
"""
@inline nondegenerate_square(b::BitIndex{Q}) where Q = iszero(UInt(b) & zero_square_bits(Q))

"""
    CliffordNumbers.sign_of_square(b::BitIndex) -> Int8

Returns the sign associated with squaring the basis blade indexed by `b` using an `Int8` as proxy:
positive signs return `Int8(1)`, negative signs return `Int8(-1)`, and zeros from degenerate
components return `Int8(0)`.
"""
sign_of_square(b::BitIndex) = (Int8(-1)^signbit_of_square(b)) * nondegenerate_square(b)

"""
    CliffordNumbers.signbit_of_mult(a::Integer, [b::Integer]) -> Bool
    CliffordNumbers.signbit_of_mult(a::BitIndex, [b::BitIndex]) -> Bool

Calculates the sign bit associated with multiplying basis elements indexed with bit indices supplied
as either integers or `BitIndex` instances. The sign bit flips when the order of `a` and `b` are
reversed, unless `a === b`. 

As with `Base.signbit()`, `true` represents a negative sign and `false` a positive sign. However,
in degenerate metrics (such as those of projective geometric algebras) the sign bit may be
irrelevant as the multiplication of those basis blades would result in zero.
"""
@inline function signbit_of_mult(a::Unsigned, b::Unsigned)
    a = abs(a) >>> 1
    sum = 0
    while !iszero(a)
        sum += count_ones(a & abs(b))
        a = a >>> 1
    end
    return isodd(sum)
end

# Account for the sign bits of signed integers
signbit_of_mult(a::Integer, b::Integer) = xor(signbit_of_mult(unsigned.(a,b)...), signbit(xor(a,b)))

@inline function signbit_of_mult(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    base_signbit = xor(signbit_of_mult(UInt(a), UInt(b)), signbit(a), signbit(b))
    return xor(base_signbit, isodious(UInt(a) & UInt(b) & negative_square_bits(Q)))
end

"""
    CliffordNumbers.nondegenerate_mult(a::T, b::T) where T<:BitIndex -> Bool

Returns `false` if the product of `a` and `b` is zero due to the squaring of a degenerate component,
`true` otherwise. This function always returns `true` if `R === 0`.
"""
@inline function nondegenerate_mult(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    return iszero(UInt(a) & UInt(b) & zero_square_bits(Q))
end

"""
    CliffordNumbers.sign_of_mult(a::T, b::T) where T<:BitIndex -> Int8

Returns an `Int8` that carries the sign associated with the multiplication of two basis blades of
Clifford/geometric algebras of the same quadratic form.
"""
function sign_of_mult(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    return Int8(-1)^signbit_of_mult(a,b) * !nondegenerate_mult(a,b)
end

#---Multiplication and duals-----------------------------------------------------------------------#
"""
    *(a::BitIndex{Q}, b::BitIndex{Q}) -> BitIndex{Q}

Returns the `BitIndex` corresponding to the basis blade resulting from the geometric product of the
basis blades indexed by `a` and `b`.
"""
@inline *(a::T, b::T) where T<:BitIndex = T(signbit_of_mult(a,b), xor(UInt(a), UInt(b)))

"""
    CliffordNumbers.has_wedge(a::BitIndex{Q}, b::BitIndex{Q}, [c::BitIndex{Q}...]) -> Bool

Returns `true` if the basis blades indexed by `a`, `b`, or any other blades `c...` have a nonzero
wedge product; `false` otherwise. This is determined by comparing all bits of the arguments (except
the sign bit) to identify any matching basis blades using bitwise AND.
"""
has_wedge(a::BitIndex{Q}, b::BitIndex{Q}) where Q = iszero((a.i << 1) & (b.i << 1))

function has_wedge(a::I, b::I, C::I...) where {Q,I<:BitIndex{Q}}
    # The sign of multiplication is irrelevant
    ab = I(xor(UInt(a), UInt(b)))
    return has_wedge(a, b) && has_wedge(ab, first(C), C[2:end]...)
end

dual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(false, typemax(UInt))
undual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(!iszero(dimension(Q) & 2), typemax(UInt))

# This calculates the "inverse" of i
# TODO: document this fully
_inv(i::T) where T<:BitIndex = T(xor(signbit(i), signbit_of_square(i)), UInt(i))

"""
    left_complement(b::BitIndex{Q}) -> BitIndex{Q}

Returns the left complement of `b`, define so that `left_complement(b) * b` generates the
pseudoscalar index of elements of the algebra `Q`.

When the left complement is applied twice, the original `BitIndex` object is returned up to a change
of sign, given by `(-1)^(grade(b) * (dimension(Q) - grade(b))). This implies that in algebras of odd 
dimension, the left complement and [right complement](@ref right_complement) are identical because
either `grade(b)` or `dimension(Q) - grade(b)` must be even. The complement is independent of the
signature of `Q`, depending only on the dimension.

Lengyel's convention for the left complement is an underbar.
"""
left_complement(b::BitIndex) = typeof(b)(false, typemax(UInt)) * _inv(b)

"""
    right_complement(b::BitIndex{Q}) -> BitIndex{Q}

Returns the right complement of `b`, define so that `b * right_complement(b)` generates the
pseudoscalar index of elements of the algebra `Q`.

When the right complement is applied twice, the original `BitIndex` object is returned up to a
change of sign, given by `(-1)^(grade(b) * (dimension(Q) - grade(b))). This implies that in algebras
of odd dimension, the [left complement](@ref left_complement) and right complement are identical
because either `grade(b)` or `dimension(Q) - grade(b)` must be even. The complement is independent
of the signature of `Q`, depending only on the dimension.

Lengyel's convention for the right complement is an overbar.
"""
right_complement(b::BitIndex)= _inv(b) * typeof(b)(false, typemax(UInt))
