module CliffordNumbersStructArraysExt

using CliffordNumbers
using StructArrays

import CliffordNumbers: AbstractCliffordNumber, KVector, Z2CliffordNumber, BladeIndex
import CliffordNumbers: BladeIndices, nblades, scalar_type, signature, grade
import CliffordNumbers.Metrics

#---Types eligible for structure-of-arrays (SoA) column storage------------------------------------#
#=
`StructArray(::Vector{KVector{...}})` would otherwise collapse to a single-column `StructVector`
holding the inner `NTuple` as one opaque field, which is useless for columnar kernels. By overriding
the `StructArrays` interface (`staticschema` / `component` / `createinstance`) for the concrete
`isbits` Clifford number types, each blade coefficient gets its own `Vector{T}` column instead.

The interface is keyed off the fully parameterized type (including `L`), so distinct grades and
algebras receive distinct, correctly-named schemas.
=#
const SoACliffordNumber = Union{KVector, Z2CliffordNumber}

#---Column names derived from the blade bit patterns-----------------------------------------------#
"""
    CliffordNumbersStructArraysExt._blade_name(Q, b::BladeIndex{Q}) -> Union{Symbol,Nothing}

Returns a `Symbol` naming the basis blade indexed by `b` (e.g. `:e12` for `e₁∧e₂`, `:scalar` for the
grade-0 element), or `nothing` when no unambiguous printable name exists (i.e. some basis vector
index falls outside `0:9`, which would make the concatenated name ambiguous).
"""
function _blade_name(Q, b::BladeIndex)
    iszero(grade(b)) && return :scalar
    bits = UInt(b)
    labels = Int[]
    for n in 0:(Metrics.dimension(Q) - 1)
        iszero(bits & (UInt(1) << n)) && continue
        push!(labels, Int(eachindex(Q)[n + 1]))
    end
    # A name like `e12` is only unambiguous when every basis index is a single nonnegative digit
    all(l -> 0 <= l <= 9, labels) || return nothing
    return Symbol(Metrics.blade_symbol(Q), join(labels))
end

"""
    CliffordNumbersStructArraysExt._soa_names(::Type{C}) -> NTuple{nblades(C),Symbol}

Computes the SoA column names for a Clifford number type `C`, ordered to match `Tuple(::C)`. Bit-
pattern names (`:scalar`, `:e12`, …) are used when they are all available and unique; otherwise the
whole schema falls back to positional names `:d1` … `:dL`.
"""
function _soa_names(::Type{C}) where C<:AbstractCliffordNumber
    Q = signature(C)
    inds = Tuple(BladeIndices(C))
    raw = map(b -> _blade_name(Q, b), inds)
    if all(!isnothing, raw) && allunique(raw)
        return map(Symbol, raw)
    end
    return ntuple(i -> Symbol(:d, i), Val(nblades(C)))
end

#---StructArrays interface-------------------------------------------------------------------------#

@generated function StructArrays.staticschema(::Type{C}) where C<:SoACliffordNumber
    NT = NamedTuple{_soa_names(C), NTuple{nblades(C), scalar_type(C)}}
    return :($NT)
end

@generated function StructArrays.component(x::C, key::Symbol) where C<:SoACliffordNumber
    names = _soa_names(C)
    ex = :(throw(ArgumentError(string($C, " has no component ", repr(key)))))
    for i in length(names):-1:1
        ex = :(key === $(QuoteNode(names[i])) ? (@inbounds Tuple(x)[$i]) : $ex)
    end
    return ex
end

StructArrays.component(x::SoACliffordNumber, i::Integer) = @inbounds Tuple(x)[i]

StructArrays.createinstance(::Type{C}, args...) where C<:SoACliffordNumber = C(args)

end
