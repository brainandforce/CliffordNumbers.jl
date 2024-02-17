"""
    CliffordNumbers.signmask(T::Type{<:Integer}, signbit::Bool = true) -> T1

Generates a signmask, or a string of bits where the only 1 bit is the sign bit. If `signbit` is set
to false, this returns a string of bits.
"""
signmask(T::Type{<:Integer}, signbit::Bool = true) = bitreverse(T(signbit))
signmask(x::Integer, signbit::Bool = true) = bitreverse(typeof(x)(signbit))

"""
    BitIndex{Q<:QuadraticForm}

A representation of an index corresponding to a basis blade of the geometric algebra with quadratic
form `Q`.
"""
struct BitIndex{Q<:QuadraticForm}
    i::UInt
    # Construct with the sign bit separately
    function BitIndex{Q}(signbit::Bool, blade::Unsigned) where Q
        x = ifelse(Q === QuadraticForm, signmask(blade), elements(Q))
        return new((blade % x) | signmask(blade, signbit))
    end
end

const GenericBitIndex = BitIndex{QuadraticForm}

#---Properties (defined separately from fields)----------------------------------------------------#

function Base.propertynames(::BitIndex; private=false)
    return private ? (:i, :signbit, :blade) : (:signbit, :blade)
end

function Base.getproperty(b::BitIndex, s::Symbol)
    s === :signbit && return !iszero(getfield(b, :i) & signmask(UInt))
    s === :blade && return getfield(b, :i) & ~signmask(UInt)
    return getfield(b, s)
end

#---Constructor tools------------------------------------------------------------------------------#
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
    return (v, !iseven(swaps))
end

function _bitindex!(Q::Type{<:QuadraticForm}, v::AbstractVector{<:Integer})
    (v, parity) = _sort_with_parity!(v)
    # Remove pairs of identical elements in v
    for n in axes(v,1)[2:end]
        # If we find a pair, set both elements to zero
        v[n] == v[n-1] && (v[[n,n-1]] .= 0)
    end
    return BitIndex{Q}(parity, sum(UInt(2)^(x-1) for x in filter!(!iszero, v); init=UInt(0)))
end

"""
    BitIndex{Q}(x::Integer...)
    BitIndex{Q}(v::AbstractVector{<:Integer})

Constructs a `BitIndex{Q}` from a list of integers that represent the basis vectors of the space.

This package uses a lexicographic convention for basis blades: in the algebra of physical space, the
basis bivectors are {e₁e₂, e₁e₃, e₂e₃}. The sign of the `BitIndex{Q}` is negative when the parity of
the basis vector permutation is odd.
"""
BitIndex{Q}(x::Integer...) where Q = _bitindex!(Q, collect(x))
BitIndex{Q}(v::AbstractVector{<:Integer}) where Q = _bitindex!(Q, deepcopy(v))

BitIndex(x) = GenericBitIndex(x)

#---Show method------------------------------------------------------------------------------------#

function Base.show(io::IO, b::BitIndex{Q}) where Q
    print(io, "-"^b.signbit, typeof(b), "(")
    found_first_vector = false
    for a in 1:min(dimension(Q), 8*sizeof(b.i) - 1)
        if !iszero(b.i & 2^(a-1))
            found_first_vector && print(io, ", ")
            print(io, a)
            found_first_vector = true
        end
    end
    print(io, ")")
end

#---Other useful functions-------------------------------------------------------------------------#

Base.signbit(b::BitIndex) = b.signbit
Base.sign(b::BitIndex) = Int8(-1)^b.signbit

Base.:-(b::BitIndex) = typeof(b)(!b.signbit, b.blade)
Base.abs(b::BitIndex) = typeof(b)(false, b.blade)

grade(b::BitIndex) = count_ones(b.blade)

#---Grade dependent sign inversion-----------------------------------------------------------------#
import Base: reverse, conj

"""
    reverse(b::BitIndex) -> BitIndex
    reverse(x::AbstractCliffordNumber{Q,T}) -> typeof(x)

Performs the reverse operation on the basis blade indexed by `b` or the Clifford number `x`. The 
sign of the reverse depends on the grade, and is positive for `g % 4 in 0:1` and negative for
`g % 4 in 2:3`.
"""
Base.reverse(b::BitIndex) = typeof(b)(xor(signbit(b), !iszero(grade(b) & 2)), b.blade)

"""
    grade_involution(b::BitIndex) -> BitIndex
    grade_involution(x::AbstractCliffordNumber{Q,T}) -> typeof(x)

Calculates the grade involution of the basis blade indexed by `b` or the Clifford number `x`. This
effectively reflects all of the basis vectors of the space along their own mirror operation, which
makes elements of odd grade flip sign.
"""
grade_involution(b::BitIndex) = typeof(b)(xor(signbit(b), !iseven(b)), b.blade)

"""
    conj(b::BitIndex) -> BitIndex
    conj(x::AbstractCliffordNumber{Q,T}) -> typeof(x)

Calculates the Clifford conjugate of the basis blade indexed by `b` or the Clifford number `x`. This
is equal to `grade_involution(reverse(x))`.
"""
Base.conj(b::BitIndex) = typeof(b)(xor(signbit(b), !iszero(grade(b)+1 & 2)), b.blade)

#---Multiplication tools---------------------------------------------------------------------------#
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
function signbit_of_mult(a::Unsigned, b::Unsigned)
    a = abs(a) >>> 1
    sum = 0
    while !iszero(a)
        sum += count_ones(a & abs(b))
        a = a >>> 1
    end
    return !iszero(sum & 1)
end

# Account for the sign bits of signed integers
signbit_of_mult(a::Integer, b::Integer) = xor(signbit_of_mult(unsigned.(a,b)...), signbit(xor(a,b)))
signbit_of_mult(a::GenericBitIndex, b::GenericBitIndex) = signbit_of_mult(a.blade, b.blade)

function signbit_of_mult(
    a::BitIndex{QuadraticForm{P,Q,R}},
    b::BitIndex{QuadraticForm{P,Q,R}}
) where {P,Q,R}
    base_signbit = xor(signbit_of_mult(a.blade, b.blade), a.signbit, b.signbit)
    iszero(Q) && return base_signbit
    # Only perform this test for pseudo-Riemannian metrics
    q = sum(UInt(2)^(n-1) for n in P .+ (1:Q); init=0)
    return xor(base_signbit, !isevil(a.blade & b.blade & q))
end

signbit_of_mult(i) = signbit_of_mult(i,i)

"""
    CliffordNumbers.sign_of_mult(a::T, b::T) where T<:BitIndex{QuadraticForm{P,Q,R}} -> Int8

Returns an `Int8` that carries the sign associated with the multiplication of two basis blades of
Clifford/geometric algebras of the same quadratic form.
"""
function sign_of_mult(
    a::BitIndex{QuadraticForm{P,Q,R}},
    b::BitIndex{QuadraticForm{P,Q,R}}
) where {P,Q,R}
    base_signbit = signbit_of_mult(a.blade, b.blade)
    # For Euclidean spaces no further processing is needed
    iszero(R) && return Int8(-1)^base_signbit
    # If any dimension squares to zero, just return zero
    r = sum(UInt(2)^(n-1) for n in (P + Q) .+ (1:R); init=0)
    return iszero(a.blade & b.blade & r) ? Int8(-1)^base_signbit : Int8(0)
end

sign_of_mult(a::GenericBitIndex, b::GenericBitIndex) = Int8(-1)^signbit_of_mult(a,b)
sign_of_mult(i) = sign_of_mult(i,i)

#---Multiplication and duals-----------------------------------------------------------------------#

function Base.:*(a::BitIndex{Q}, b::BitIndex{Q}) where Q
    return BitIndex{Q}(signbit_of_mult(a,b), xor(a.blade, b.blade))
end

dual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(false, typemax(UInt))
undual(b::BitIndex{Q}) where Q = b * BitIndex{Q}(!iszero(dimension(Q) & 2), typemax(UInt))

#---Clifford number iteration----------------------------------------------------------------------#
"""
    BitIndices{Q<:QuadraticForm,C<:AbstractCliffordNumber{Q,<:Any}} <: AbstractVector{BitIndex{Q}}

Represents a range of valid `BitIndex` objects for the nonzero components of a given multivector 
with quadratic form `Q`.

For a generic `AbstractCliffordNumber{Q}`, this returns `BitIndices{CliffordNumber{Q}}`, which
contains all possible indices for a multivector associated with the quadratic form `Q`. This may 
also be constructed with `BitIndices(Q)`.

For sparse representations, such as `KVector{K,Q}`, the object only contains the indices of the
nonzero elements of the multivector.

# Construction

`BitIndices` can be constructed by calling the type constructor with either the multivector or its
type.

# Indexing

`BitIndices` always uses one-based indexing like most Julia arrays. Although it is more natural in
the dense case to use zero-based indexing, as the basis blades are naturally encoded in the indices
for the dense representation of `CliffordNumber`, one-based indexing is used by the tuples which
contain the data associated with this package's implementations of Clifford numbers.

# Interfaces for new subtypes of `AbstractCliffordNumber`

When defining the behavior of `BitIndices` for new subtypes `T` of `AbstractCliffordNumber`, 
`Base.getindex(::BitIndices{Q,T}, i::Integer)` should be defined so that all indices of T that are
not constrained to be zero are returned.
"""
struct BitIndices{Q,C<:AbstractCliffordNumber{Q,<:Any}} <: AbstractVector{BitIndex{Q}}
end

BitIndices{Q}(::Type{C}) where {Q,C<:AbstractCliffordNumber} = BitIndices{Q,C}()
BitIndices(C::Type{<:AbstractCliffordNumber{Q,<:Any}}) where Q = BitIndices{Q,C}()

BitIndices(x::AbstractCliffordNumber{Q,<:Any}) where Q = BitIndices{Q,typeof(x)}()

Base.size(::BitIndices{Q,C}) where {Q,C} = tuple(length(C))

function Base.getindex(b::BitIndices{Q}, i::Integer) where Q
    @boundscheck checkbounds(b, i)
    return BitIndex{Q}(signbit(i-1), unsigned(i-1))
end

#---Range of valid indices for CliffordNumber------------------------------------------------------#

Base.keys(x::AbstractCliffordNumber) = keys(typeof(x))  # only need to define on types
Base.keys(::Type{T}) where T<:AbstractCliffordNumber{<:Any,<:Any} = BitIndices(T)

#---Indexing an AbstractCliffordNumber with BitIndices of a type-----------------------------------#

function Base.getindex(x::AbstractCliffordNumber{Q}, b::BitIndices{Q,T}) where {Q,T}
    data = ntuple(i -> x[b[i]], Val(length(T)))
    return T(data)
end
