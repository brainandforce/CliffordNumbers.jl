module CliffordNumbersChainRulesCoreExt

using CliffordNumbers

import ChainRulesCore: NoTangent, ProjectTo, frule, rrule

# Generic rrule for all AbstractCliffordNumber instances
# TODO: check that this works for a variety of cases
# It does not seem to get gradients correct for f(x) = x^2
# But it does seem to work for f(x) = x*x
function rrule(::Type{T}, t::Tuple) where T<:AbstractCliffordNumber
    CliffordNumber_pullback(∇x) = (NoTangent(), Tuple(∇x))
    return (T(t), CliffordNumber_pullback)
end

end
