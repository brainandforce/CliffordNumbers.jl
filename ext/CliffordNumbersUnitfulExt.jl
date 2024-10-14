module CliffordNumbersUnitfulExt

using CliffordNumbers
using Unitful

import CliffordNumbers: ∧, ∨, ⨼, ⨽, dot, ×, ⨰
import Unitful: AbstractQuantity, AffineError, AffineQuantity

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
