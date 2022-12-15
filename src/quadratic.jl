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
    APS

The algebra of physical space, Cl(3,0,0).
An alias for `QuadraticForm{3,0,0}`.
"""
const APS = QF{3,0,0}

"""
    STAplus

Spacetime algebra with a mostly positive signature (mathematician's convention), Cl(3,1,0).
An alias for `QuadraticForm{3,1,0}`.
"""
const STAplus = QF{3,1,0}

"""
    STAminus

Spacetime algebra with a mostly negative signature (particle physicist's convention), Cl(1,3,0).
An alias for `QuadraticForm{1,3,0}`.
"""
const STAminus = QF{1,3,0}

"""
    PGA3D

Projective geometric algebra of 2 dimensions and 1 projective dimension, Cl(2,0,1).
An alias for `QuadraticForm{2,0,1}`.
"""
const PGA2D = QF{2,0,1}

"""
    PGA3D

Projective geometric algebra of 3 dimensions and 1 projective dimension, Cl(3,0,1).
An alias for `QuadraticForm{3,0,1}`.
"""
const PGA3D = QF{3,0,1}
