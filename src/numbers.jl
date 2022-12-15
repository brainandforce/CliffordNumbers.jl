"""
    CliffordNumber{Cl,T,L}

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Cl`, 
element type `T`, and length `L` (which depends entirely on `Cl`).

The coefficients are ordered by taking advantage of the natural binary structure of the basis. The
grade of an element is given by the Hamming weight of its index. For the algebra of physical space,
the order is: 1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃ = i. This order allows for more aggressive SIMD
optimization when calculating the geometric product.
"""
struct CliffordNumber{Cl<:QuadraticForm,T<:Number,L} <: Number
    data::NTuple{L,T}
    function CliffordNumber{Cl,T,L}(x) where {Cl,T,L}
        sz = elements(Cl)
        @assert length(x) == L == sz string(
            "Incorrect number of components: multivectors of ", Cl, " have ", sz, " components."
        )
        return new{Cl,T,L}(x)
    end
end

#---Constructors----------------------------------------------------------------------------------#

CliffordNumber{Cl,T}(x::NTuple{L,<:Number}) where {Cl,T,L} = CliffordNumber{Cl,T,L}(x)
CliffordNumber{Cl,T}(x::Vararg{<:Number,L}) where {Cl,T,L} = CliffordNumber{Cl,T,L}(x)

# A constructor similar to `ntuple(::Function)`
# However, it deals with the offset indexing
function CliffordNumber{Cl,T}(f::Function) where {Cl,T}
    L = elements(Cl)
    return CliffordNumber{Cl,T,L}(ntuple(i -> f(i-1), Val{L}()))
end

# Promote to a common type first 

function CliffordNumber{Cl}(x::NTuple{L,<:Number}) where {Cl,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Cl,T,L}(promote(x...))
end

function CliffordNumber{Cl}(x::Vararg{<:Number,L}) where {Cl,L}
    T = promote_type(typeof.(x)...)
    return CliffordNumber{Cl,T,L}(promote(x...))
end

function CliffordNumber{Cl}(f::Function) where {Cl}
    L = elements(Cl)
    data = ntuple(i -> f(i-1), Val{L}())
    return CliffordNumber{Cl,eltype(data),L}(data)
end

function CliffordNumber{Cl,T}(x::Real) where {Cl,T}
    L = elements(Cl)
    return CliffordNumber{Cl,T,L}(ntuple(i -> T(isone(i) * x), Val{L}()))
end

function CliffordNumber{Cl,T1}(x::Complex{T2}) where {Cl,T1<:Real,T2<:Real}
    L = elements(Cl)
    T = promote_type(T1,T2)
    data = ntuple(Val{L}()) do i
        i == 1 && return T(real(x))
        i == L && return T(imag(x))
    end
    return CliffordNumber{Cl,T,L}(data)
end

CliffordNumber{Cl}(x::Number) where Cl = CliffordNumber{Cl,typeof(x)}(x)

#---Number of elements----------------------------------------------------------------------------#

import Base.length
# This is equal to the `L` parameter
length(::Type{<:CliffordNumber{Cl}}) where Cl = elements(Cl)
length(m::CliffordNumber{Cl}) where Cl = length(typeof(m))

#---Get type parameters---------------------------------------------------------------------------#

Base.eltype(::Type{<:CliffordNumber{Cl,T}}) where {Cl,T} = T
algebra(::Type{<:CliffordNumber{Cl}}) where Cl = Cl
algebra(::CliffordNumber{Cl}) where Cl = Cl

#---Generate zero and identity elements-----------------------------------------------------------#

import Base: zero, one, oneunit

zero(S::Type{<:CliffordNumber{Cl,T}}) where {Cl,T} = S(i -> zero(T))
oneunit(S::Type{<:CliffordNumber{Cl,T}}) where {Cl,T} = S(iszero)
one(S::Type{<:CliffordNumber{Cl,T}}) where {Cl,T} = oneunit(S)

function pseudoscalar(::Type{<:CliffordNumber{Cl,T}}) where {Cl,T}
    L = elements(Cl)
    return CliffordNumber{Cl,T,L}(ntuple(isequal(L), Val{L}()))
end

#---Constructors using just the quadratic forms---------------------------------------------------#

for fn in (:zero, :one, :oneunit, :pseudoscalar)
    # Default to Bools since they are promoted to any wider type
    @eval $fn(::Type{CliffordNumber{Cl}}) where Cl = $fn(CliffordNumber{Cl,Bool})
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

summary(io::IO, m::CliffordNumber) = println(io, typeof(m), ":")

#=
function show(io::IO, ::MIME"text/plain", m::CliffordNumber{Cl}) where Cl
    summary(io, m)
    # Loop through all the grades
    for n in grades(Cl)
        # Find all numbers with specific Hamming weights
        inds = findall(x -> hamming_weight(x) == n, grades(Cl))
        for i in inds
            if !iszero(m[i])
                print(io, sign(m[i]) > 0 ? "+" : "-", " ")
                print(io, abs(m[i]), to_basis_str(Cl, i), " ")
            end
        end
    end
end
=#
