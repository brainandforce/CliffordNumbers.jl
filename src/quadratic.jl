"""
    CliffordNumbers.QuadratricForm

Represents a quadratric with `P` dimensions which square to +1, `Q` dimensions which square to -1,
and `R` dimensions which square to 0, in that order.

By convention, this type is used as a tag, and is never instantiated.
"""
struct QuadraticForm{P,Q,R}
end

dimension(::Type{QuadraticForm{P,Q,R}}) where {P,Q,R} = (P + Q + R)
isdegenerate(::Type{QuadraticForm{P,Q,R}}) where {P,Q,R} = !iszero(R)
iseuclidean(::Type{QuadraticForm{P,Q,R}}) where {P,Q,R} = (iszero(Q) && iszero(R))

elements(Q::Type{<:QuadraticForm}) = 2^dimension(Q)
grades(Q::Type{<:QuadraticForm}) = 0:dimension(Q)

"""
    sign(::Type{QuadraticForm{P,Q,R}}, i::Integer) -> Int8

Gets the sign associated with dimension `i` of a quadratric form.
"""
function sign(::Type{QuadraticForm{P,Q,R}}, i::Integer) where {P,Q,R} 
    return Int8(-1)^(i in P .+ (1:Q)) * !(i in (P + Q) .+ (1:R))
end

#---Special geometric algebras---------------------------------------------------------------------#
"""
    CliffordNumbers.QFComplex

The quadratic form with zero dimensions, `QuadraticForm{0,0,0}`, isomorphic to the real numbers.
"""
const QFReal = QuadraticForm{0,0,0}

"""
    CliffordNumbers.QFComplex

The quadratic form with one dimension squaring to -1, `QuadraticForm{0,1,0}`. This generates a
Clifford algebra isomorphic to the complex numbers.
"""
const QFComplex = QuadraticForm{0,1,0}

"""
    CliffordNumbers.QFPositiveDefinite{P} (alias for QuadraticForm{P,0,0})

A positive-definite quadratic form with `P` dimensions.
"""
const QFPositiveDefinite{P} = QuadraticForm{P,0,0}

"""
    CliffordNumbers.QFNondegenerate{P,Q} (alias for QuadraticForm{P,Q,0})

Represents a non-degenerate quadratic form (one without dimensions which square to 0).
"""
const QFNondegenerate{P,Q} = QuadraticForm{P,Q,0}

"""
    CliffordNumbers.QFExterior{R} (alias for QuadraticForm{0,0,R})

The quadratic form associated with an exterior algebra, which can be thought of as a Clifford 
algebra where all dimensions are degenerate.
"""
const QFExterior{R} = QuadraticForm{0,0,R}

"""
    APS

The algebra of physical space, Cl(3,0,0). An alias for `QuadraticForm{3,0,0}`.
"""
const APS = QuadraticForm{3,0,0}

"""
    STA

Spacetime algebra with a mostly negative signature (particle physicist's convention), Cl(1,3,0). An
alias for `QuadraticForm{1,3,0}`.

The negative signature is used by default to distinguish this algebra from conformal geometric
algebras, which use a mostly positive signature by convention.
"""
const STA = QuadraticForm{1,3,0}

# These are functions because you can't perform the necessary addition in CGAs with a const.

"""
    VGA(D) -> Type{QuadraticForm{D,0,0}}

Creates the type of a quadratic form associated with a vector/vanilla geometric algebra (VGA) of
dimension `D`.

For reasons of type stability, avoid calling this function without constant arguments.
"""
VGA(D) = QuadraticForm{D,0,0}

"""
    PGA(D) -> Type{QuadraticForm{D,0,1}}

Creates the type of a quadratic form associated with a projective geometric algebra (PGA) of
dimension `D`.

For reasons of type stability, avoid calling this function without constant arguments.
"""
PGA(D) = QuadraticForm{D,0,1}

"""
    CGA(D) -> Type{QuadraticForm{D+1,1,0}}

Creates the type of a quadratic form associated with a conformal geometric algebra (CGA) of
dimension `D`.

For reasons of type stability, avoid calling this function without constant arguments.
"""
CGA(D) = QuadraticForm{D+1,1,0}

#---Show methods for special quadratic forms-------------------------------------------------------#

# Module prefix is left here just to avoid some incorrect warnings with VS Code extension
Base.show(io::IO, ::Type{VGA(D)}) where D = print(io, VGA, "($D)")
Base.show(io::IO, ::Type{PGA(D)}) where D = print(io, PGA, "($D)")
Base.show(io::IO, ::Type{QuadraticForm{D,1,0}}) where D = print(io, CGA, "($(D-1))")

# Functionality stolen from Julia base/show.jl around line 927
function Base.show(io::IO, ::MIME"text/plain", T::Type{VGA(D)}) where D
    show(IOContext(io, :compact => true), T)
    if !(get(io, :compact, false)::Bool)
        printstyled(io, " (alias for $QuadraticForm{$D,0,0})"; color = :light_black)
    end
end

function Base.show(io::IO, ::MIME"text/plain", T::Type{PGA(D)}) where D
    show(IOContext(io, :compact => true), T)
    if !(get(io, :compact, false)::Bool)
        printstyled(io, " (alias for $QuadraticForm{$D,0,1})"; color = :light_black)
    end
end

function Base.show(io::IO, ::MIME"text/plain", T::Type{QuadraticForm{D,1,0}}) where D
    show(IOContext(io, :compact => true), T)
    if !(get(io, :compact, false)::Bool)
        printstyled(io, " (alias for $QuadraticForm{$D,0,0})"; color = :light_black)
    end
end
