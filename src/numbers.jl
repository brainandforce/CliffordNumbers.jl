#---Dense representation of Clifford numbers-------------------------------------------------------#
"""
    CliffordNumber{Q,T,L} <: AbstractCliffordNumber{Q,T}

A dense multivector (or Clifford number), with quadratic form `Q`, element type `T`, and length `L`
(which depends entirely on `Q`).

The coefficients are ordered by taking advantage of the natural binary structure of the basis. The
grade of an element is given by the Hamming weight of its index. For the algebra of physical space,
the order is: 1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃ = i. This order allows for more aggressive SIMD
optimization when calculating the geometric product.
"""
struct CliffordNumber{Q<:QuadraticForm,T<:BaseNumber,L} <: AbstractCliffordNumber{Q,T}
    data::NTuple{L,T}
    function CliffordNumber{Q,T,L}(x) where {Q,T,L}
        sz = elements(Q)
        @assert length(x) == L == sz string(
            "Incorrect number of components: multivectors of ", Q, " have ", sz, " components."
        )
        return new{Q,T,L}(x)
    end
end

#---Constructors-----------------------------------------------------------------------------------#

CliffordNumber{Q,T}(x::NTuple{L,<:BaseNumber}) where {Q,T,L} = CliffordNumber{Q,T,L}(x)
CliffordNumber{Q,T}(x::Vararg{BaseNumber,L}) where {Q,T,L} = CliffordNumber{Q,T,L}(x)

# Constructors similar to `ntuple(::Function)`
# However, it deals with the offset indexing
function CliffordNumber{Q,T,L}(f::Function) where {Q,T,L}
    return CliffordNumber{Q,T,L}(ntuple(i -> f(i-1), Val{L}()))
end

CliffordNumber{Q,T}(f::Function) where {Q,T} = CliffordNumber{Q,T,elements(Q)}(f)

# Promote to a common type first 

function CliffordNumber{Q}(x::NTuple{L,<:BaseNumber}) where {Q,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Q,T,L}(promote(x...))
end

function CliffordNumber{Q}(x::Vararg{BaseNumber,L}) where {Q,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Q,T,L}(promote(x...))
end

function CliffordNumber{Q}(f::Function) where {Q}
    L = elements(Q)
    data = ntuple(i -> f(i-1), Val{L}())
    return CliffordNumber{Q,eltype(data),L}(data)
end

function CliffordNumber{Q,T,L}(x::Real) where {Q,T,L}
    return CliffordNumber{Q,T,L}(ntuple(i -> T(isone(i) * x), Val{L}()))
end

CliffordNumber{Q,T}(x::Real) where {Q,T} = CliffordNumber{Q,T,elements(Q)}(x)

function CliffordNumber{Q,T1}(x::Complex{T2}) where {Q,T1<:Real,T2<:Real}
    L = elements(Q)
    T = promote_type(T1,T2)
    data = ntuple(Val{L}()) do i
        i == 1 && return T(real(x))
        i == L && return T(imag(x))
    end
    return CliffordNumber{Q,T,L}(data)
end

CliffordNumber{Q}(x::BaseNumber) where Q = CliffordNumber{Q,typeof(x)}(x)

#---Number of elements-----------------------------------------------------------------------------#

import Base.length
# This is equal to the `L` parameter
length(::Type{<:CliffordNumber{Q}}) where Q = elements(Q)
length(m::CliffordNumber) = length(typeof(m))

nonzero_grades(::Type{<:CliffordNumber{Q}}) where Q = 0:dimension(Q)

#---Default BitIndices construction should include all possible BitIndex objects-------------------#

BitIndices{Q}() where Q = BitIndices{Q,CliffordNumber{Q}}()
BitIndices(Q::Type{<:QuadraticForm}) = BitIndices{Q}()

#---Clifford number indexing-----------------------------------------------------------------------#

Base.getindex(x::CliffordNumber{Q}, b::BitIndex{Q}) where Q = sign(b) * x.data[b.blade + 1]
Base.getindex(x::CliffordNumber, b::GenericBitIndex) = sign(b) * x.data[b.blade % length(x) + 1]

#---Generate zero and identity elements------------------------------------------------------------#

import Base: zero, one, oneunit

zero(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S(i -> zero(T))
oneunit(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S(iszero)
one(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = oneunit(S)

function pseudoscalar(::Type{<:CliffordNumber{Q,T}}) where {Q,T}
    L = elements(Q)
    return CliffordNumber{Q,T,L}(ntuple(isequal(L), Val{L}()))
end

pseudoscalar(m::CliffordNumber) = pseudoscalar(typeof(m))

"""
    isscalar(m::CliffordNumber)

Determines whether the Clifford number `m` is a scalar, meaning that it has no components with
grades above zero.
"""
isscalar(m::CliffordNumber) = all(iszero, m.data[2:end])
isscalar(m::CliffordNumber{QuadraticForm{0,0,0}}) = true

"""
    ispseudoscalar(m::CliffordNumber)

Determines whether the Clifford number `m` is a pseudoscalar, meaning that it has no components with
grades below the one equal to the dimension of the space.
"""
ispseudoscalar(m::CliffordNumber) = all(iszero, m.data[1:end-1])

#---Constructors using just the quadratic forms----------------------------------------------------#

for fn in (:zero, :one, :oneunit, :pseudoscalar)
    # Default to Bools since they are promoted to any wider type
    @eval $fn(::Type{CliffordNumber{Q}}) where Q = $fn(CliffordNumber{Q,Bool})
    @eval $fn(Q::Type{<:QuadraticForm}) = $fn(CliffordNumber{Q,Bool})
end

#---Show methods-----------------------------------------------------------------------------------#
import Base: show, summary

function show(io::IO, m::CliffordNumber)
    print(io, "CliffordNumber{", QuadraticForm(m), ",", eltype(m), "}(", join(m.data, ", "), ")")
end

function to_basis_str(
    Cl::Type{QuadraticForm{P,Q,R}},
    i::Integer;
    pseudoscalar=""
) where {P,Q,R}
    if i >= elements(Cl)
        error("i = $i exceeds maximum number of elements for $Cl")
    elseif i == elements(Cl) - 1 && !isempty(pseudoscalar)
        return string(pseudoscalar)
    end
    return join(
        [("e" * subscript_string(n+1))^!iszero(typeof(i)(2)^n & i) for n in 0:(dimension(Cl) - 1)]
    )
end

function summary(io::IO, x::CliffordNumber)
    println(io, "CliffordNumber{", QuadraticForm(x), ",", eltype(x), "}:")
end

function show(io::IO, ::MIME"text/plain", x::CliffordNumber{Q}) where Q
    summary(io, x)
    b = BitIndices(x)
    # Flag to mark when we've *found the first nonzero* element
    ffn = false
    # Print the scalar component first
    if !iszero(x[b[1]])
        print(io, x[b[1]])
        ffn = true
    end
    # Loop through all the grades
    for n in 1:dimension(Q)
        # Find all numbers with specific Hamming weights
        inds = findall(x -> count_ones(x) == n, 0:(elements(Q) - 1))
        for i in inds
            if !iszero(x[b[i]])
                print(
                    io, " "^ffn, (sign(x[b[i]]) > 0 ? "+"^ffn : "-"), " "^ffn, abs(x[b[i]]), 
                    "*"^(x[b[i]] isa Bool), to_basis_str(Q, i-1)
                )
                ffn = true
            end
        end
    end
    # If we got through the whole process, just print the zero
    !ffn && print(io, x[b[1]])
end
