#---Abstract type for all Clifford numbers---------------------------------------------------------#
"""
    AbstractCliffordNumber{Q,T} <: Number

An element of a Clifford algebra, often referred to as a multivector, with quadratic form `Q` and
element type `T`.
"""
abstract type AbstractCliffordNumber{Q<:QuadraticForm,T<:BaseNumber} <: Number
end

#---Get type parameters----------------------------------------------------------------------------#

Base.eltype(::Type{<:AbstractCliffordNumber{Q,T}}) where {Q,T} = T
algebra(::Type{<:AbstractCliffordNumber{Q}}) where Q = Q
algebra(::AbstractCliffordNumber{Q}) where Q = Q
