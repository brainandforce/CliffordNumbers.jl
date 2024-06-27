module CliffordNumbersStaticArraysCoreExt

using CliffordNumbers
using StaticArraysCore

function StaticArraysCore.similar_type(::Type{C}, args...) where C<:AbstractCliffordNumber
    return CliffordNumbers.similar_type(C, args...)
end

function StaticArraysCore.similar_type(x::AbstractCliffordNumber, args...) 
    return CliffordNumbers.similar_type(x, args...)
end

end
