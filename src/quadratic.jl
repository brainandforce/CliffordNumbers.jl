"""
    CliffordNumbers.QuadratricForm

Represents a quadratric with `P` dimensions which square to +1, `Q` dimensions which square to -1,
and `R` dimensions which square to 0, in that order.

By convention, this type is used as a tag, and is never instantiated.
"""
struct QuadraticForm{P,Q,R}
end

const QF = QuadraticForm

dimension(::Type{QuadraticForm{P,Q,R}}) where {P,Q,R} = (P + Q + R)
elements(Q::Type{<:QuadraticForm}) = 2^dimension(Q)
grades(Q::Type{<:QuadraticForm}) = 0:dimension(Q)

"""
    sign(::Type{QuadraticForm{P,Q,R}}, i::Integer) -> Int8

Gets the sign associated with dimension `i` of a quadratric form.
"""
function Base.sign(::Type{QuadraticForm{P,Q,R}}, i::Integer) where {P,Q,R} 
    return Int8(-1)^(i in P .+ (1:Q)) * !(i in (P + Q) .+ (1:R))
end

#---Special geometric algebras---------------------------------------------------------------------#

"""
    VGA{D} (alias for `QuadraticForm{D,0,0}`)

Alias for vector/vanilla geometric algebras, which have positive-definite signatures.
"""
const VGA{D} = QuadraticForm{D,0,0}

"""
    PGA{D} (alias for `QuadraticForm{D,0,1}`)

Alias for projective geometric algebras, which have positive-definite signatures except in one
dimension, which is degenerate (squares to zero).
"""
const PGA{D} = QuadraticForm{D,0,1}

"""
    APS

The algebra of physical space, Cl(3,0,0).
An alias for `QuadraticForm{3,0,0}`.
"""
const APS = QF{3,0,0}

"""
    STA

Spacetime algebra with a mostly negative signature (particle physicist's convention), Cl(1,3,0).
An alias for `QuadraticForm{1,3,0}`.

The negative signature is used by default to distinguish this algebra from conformal geometric
algebras, which use a mostly positive signature by convention.
"""
const STAminus = QF{1,3,0}
