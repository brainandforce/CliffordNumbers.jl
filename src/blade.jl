"""
    Blade{B,Q,T} <: AbstractCliffordNumber{Q,T}

Represents a scaled basis blade of the algebra indexed by `B`, a `BitIndex` object.
"""
struct Blade{B,Q,T<:BaseNumber} <: AbstractCliffordNumber{Q,T}
    data::T
    function Blade{B,Q,T}(x)
        @assert B isa BitIndex{Q} "B parameter must be a BitIndex{$Q} instance (got $B)."
        return new(x)
    end
end

Blade{B,Q}(x::T) = Blade{B,Q,T}(x)
Blade{B}(x::T) = Blade{B,signature(B),T}(x)

nblades(::Type{<:Blade}) = 1
grade(::Type{<:Blade{B}}) where B = grade(B)
nonzero_grades(::Type{<:Blade{B}}) where B = grade(B):1:grade(B)

@inline getindex(x::Blade{B,Q}, i::BitIndex{Q}) where {B,Q} = ifelse(i === B, x, zero(x))
