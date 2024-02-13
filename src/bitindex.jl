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

#---Clifford number indexing-----------------------------------------------------------------------#

Base.getindex(x::CliffordNumber{Q}, b::BitIndex{Q}) where Q = sign(b) * x.data[b.blade + 1]
Base.getindex(x::CliffordNumber, b::GenericBitIndex) = sign(b) * x.data[b.blade % length(x) + 1]

function Base.getindex(m::CliffordNumber, i::Integer)
    # Throw the correct BoundsError for out of bounds
    return i in 0:length(m)-1 ? m.data[i+1] : throw(BoundsError(m, i))
end

#---Clifford number iteration----------------------------------------------------------------------#
"""
    BitIndices{Q<:QuadraticForm}

Represents a range of valid `BitIndex` objects for a given quadratic form.

For sparse representations, such as `KVector{K...}`, this is not the most efficient way to iterate
through all elements, as it includes indices that are known to be zero.
"""
struct BitIndices{Q<:QuadraticForm}
end

"""
    BitIndices(x::AbstractCliffordNumber{Q}) -> BitIndices{Q}()
    BitIndices(T::Type{<:AbstractCliffordNumber{Q}}) -> BitIndices{Q}()

Constructs a `BitIndices` object associated with a Clifford number or its type.
"""
BitIndices(::Type{<:AbstractCliffordNumber{Q}}) where Q = BitIndices{Q}()
BitIndices(::AbstractCliffordNumber{Q}) where Q = BitIndices{Q}()
BitIndices(::Type{Q}) where Q = BitIndices{Q}()

Base.length(::BitIndices{Q}) where Q = elements(Q)
Base.first(::BitIndices{Q}) where Q = BitIndex{Q}(false, UInt(0))
Base.last(::BitIndices{Q}) where Q = BitIndex{Q}(false, typemax(UInt) % elements(Q))
Base.eltype(::Type{BitIndices{Q}}) where Q = BitIndex{Q}

function Base.iterate(::BitIndices{Q}, i::Integer = 0) where Q
    return 0 <= i < elements(Q) ? (BitIndex{Q}(false, UInt(i)), i+1) : nothing
end

Base.getindex(::BitIndices{Q}, i::Integer) where Q = BitIndex{Q}(signbit(i), unsigned(i))

#---Range of valid indices for CliffordNumber------------------------------------------------------#

Base.keys(x::AbstractCliffordNumber) = keys(typeof(x))  # only need to define on types
Base.keys(::Type{<:CliffordNumber{Q}}) where Q = BitIndices{Q}()
