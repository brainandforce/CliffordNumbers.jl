module CliffordNumbersUnitfulExt

using CliffordNumbers
using Unitful

import CliffordNumbers: ∧, ∨, ⨼, ⨽, dot, ×, ⨰
import Unitful: AbstractQuantity, AffineError, AffineQuantity

#---Fixing printing methods------------------------------------------------------------------------#

# Probably not needed, but included anyway
Unitful.BracketStyle(::Type{<:AbstractCliffordNumber}) = Unitful.RoundBrackets()

function Base.print(io::IO, x::AbstractQuantity{<:AbstractCliffordNumber})
    # Functions with qualified names are part of the internal Unitful API.
    # This should be fine for now, but may break unexpectedly.
    if unit(x) isa Unitful.Units{()}
        print(io, x.val, " (no units)")
    else
        print(io, '(', x.val, ')')
        # This one is based on SI convention and should not change. I hope.
        Unitful.has_unit_spacing(unit(x)) && print(io, ' ')
        print(io, unit(x))
    end
end

Base.show(io::IO, ::MIME"text/plain", x::Quantity{<:AbstractCliffordNumber}) = print(io, x)

#---Definitions for different products-------------------------------------------------------------#
∧(x::AbstractQuantity, y::AbstractQuantity) = Quantity(x.val ∧ y.val, unit(x) * unit(y))

function ∧(x::Number, y::AbstractQuantity)
    y isa AffineQuantity &&
        throw(AffineError("an invalid operation was attempted with affine quantities: $x ∧ $y"))
    return Quantity(x ∧ y.val, unit(y))
end

function ∧(x::AbstractQuantity, y::Number)
    x isa AffineQuantity &&
        throw(AffineError("an invalid operation was attempted with affine quantities: $x ∧ $y"))
    return Quantity(x.val ∧ y, unit(x))
end

end
