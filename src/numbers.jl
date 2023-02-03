"""
    CliffordNumber{Q,T,L}

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Cl`, 
element type `T`, and length `L` (which depends entirely on `Cl`).

The coefficients are ordered by taking advantage of the natural binary structure of the basis. The
grade of an element is given by the Hamming weight of its index. For the algebra of physical space,
the order is: 1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃ = i. This order allows for more aggressive SIMD
optimization when calculating the geometric product.
"""
struct CliffordNumber{Q<:QuadraticForm,T<:Union{Real,Complex},L} <: Number
    data::NTuple{L,T}
    function CliffordNumber{Q,T,L}(x) where {Q,T,L}
        sz = elements(Cl)
        @assert length(x) == L == sz string(
            "Incorrect number of components: multivectors of ", Q, " have ", sz, " components."
        )
        return new{Q,T,L}(x)
    end
end

#---Constructors----------------------------------------------------------------------------------#

CliffordNumber{Q,T}(x::NTuple{L,<:Number}) where {Q,T,L} = CliffordNumber{Q,T,L}(x)
CliffordNumber{Q,T}(x::Vararg{<:Number,L}) where {Q,T,L} = CliffordNumber{Q,T,L}(x)

# Constructors similar to `ntuple(::Function)`
# However, it deals with the offset indexing
function CliffordNumber{Q,T,L}(f::Function) where {Q,T,L}
    return CliffordNumber{Q,T,L}(ntuple(i -> f(i-1), Val{L}()))
end

CliffordNumber{Q,T}(f::Function) where {Q,T} = CliffordNumber{Q,T,elements(Cl)}(f)

# Promote to a common type first 

function CliffordNumber{Q}(x::NTuple{L,<:Number}) where {Q,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Q,T,L}(promote(x...))
end

function CliffordNumber{Q}(x::Vararg{<:Number,L}) where {Q,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Q,T,L}(promote(x...))
end

function CliffordNumber{Q}(f::Function) where {Q}
    L = elements(Cl)
    data = ntuple(i -> f(i-1), Val{L}())
    return CliffordNumber{Cl,eltype(data),L}(data)
end

function CliffordNumber{Q,T,L}(x::Real) where {Q,T,L}
    return CliffordNumber{Q,T,L}(ntuple(i -> T(isone(i) * x), Val{L}()))
end

CliffordNumber{Q,T}(x::Real) where {Q,T} = CliffordNumber{Q,T,elements(Q)}(x)

function CliffordNumber{Q,T1}(x::Complex{T2}) where {Q,T1<:Real,T2<:Real}
    L = elements(Cl)
    T = promote_type(T1,T2)
    data = ntuple(Val{L}()) do i
        i == 1 && return T(real(x))
        i == L && return T(imag(x))
    end
    return CliffordNumber{Q,T,L}(data)
end

CliffordNumber{Q}(x::Number) where Q = CliffordNumber{Q,typeof(x)}(x)

#---Number of elements----------------------------------------------------------------------------#

import Base.length
# This is equal to the `L` parameter
length(::Type{<:CliffordNumber{Q}}) where Q = elements(Q)
length(m::CliffordNumber{Q}) where Q = length(typeof(m))

#---Get type parameters---------------------------------------------------------------------------#

Base.eltype(::Type{<:CliffordNumber{Q,T}}) where {Q,T} = T
algebra(::Type{<:CliffordNumber{Q}}) where Q = Q
algebra(::CliffordNumber{Q}) where Q = Q

#---Generate zero and identity elements-----------------------------------------------------------#

import Base: zero, one, oneunit

zero(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S(i -> zero(T))
oneunit(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = S(iszero)
one(S::Type{<:CliffordNumber{Q,T}}) where {Q,T} = oneunit(S)

function pseudoscalar(::Type{<:CliffordNumber{Q,T}}) where {Q,T}
    L = elements(Cl)
    return CliffordNumber{Q,T,L}(ntuple(isequal(L), Val{L}()))
end

pseudoscalar(m::CliffordNumber) = pseudoscalar(typeof(m))

#---Constructors using just the quadratic forms---------------------------------------------------#

for fn in (:zero, :one, :oneunit, :pseudoscalar)
    # Default to Bools since they are promoted to any wider type
    @eval $fn(::Type{CliffordNumber{Q}}) where Q = $fn(CliffordNumber{Cl,Bool})
    @eval $fn(Cl::Type{<:QuadraticForm}) = $fn(CliffordNumber{Cl,Bool})
end

#---Show methods----------------------------------------------------------------------------------#
import Base: show, summary

function show(io::IO, m::CliffordNumber)
    print(io, "CliffordNumber{", algebra(m), ",", eltype(m), "}(", join(m.data, ", "), ")")
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

function summary(io::IO, m::CliffordNumber)
    println(io, "CliffordNumber{", algebra(m), ",", eltype(m), "}:")
end

function show(io::IO, ::MIME"text/plain", m::CliffordNumber{Q}) where Q
    summary(io, m)
    # For a zero multivector, just print zero
    # Make sure to cover the signed zero case for floating point elements
    iszero(m) && print(io, m[0])
    # Flag to mark when we've *found the first nonzero* element
    ffn = false
    # Print the scalar component first
    if !iszero(m[0])
        print(io, m[0])
        ffn = true
    end
    # Loop through all the grades
    for n in 1:dimension(Cl)
        # Find all numbers with specific Hamming weights
        inds = findall(x -> hamming_weight(x) == n, 0:(elements(Cl) - 1))
        for i in inds .- 1
            if !iszero(m[i])
                print(io, " "^ffn, sign(m[i]) > 0 ? "+"^ffn : "-")
                print(io, " "^ffn, abs(m[i]), to_basis_str(Cl, i, pseudoscalar="i"))
                ffn = true
            end
        end
    end
end
