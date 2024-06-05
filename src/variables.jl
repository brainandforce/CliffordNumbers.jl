#---Generate basis elements for convenience--------------------------------------------------------#

basis_vars_names(Q, prefix) = [Symbol(prefix, x) for x in eachindex(Q)]

function basis_vars_names(Q::CGA, prefix = Metrics.blade_symbol(Q))
    first_names = [Symbol(prefix, :_n), Symbol(prefix, :_p)]
    return append!(first_names, [Symbol(prefix, x) for x in 1:dimension(Q) - 2])
end

function basis_vars_expr(
    ::Val{Q},
    ::Type{T} = Int,
    prefix = Metrics.blade_symbol(Q)
) where {Q,T<:BaseNumber}
    vars = basis_vars_names(Q, prefix)
    exs = map(enumerate(vars)) do (n,var)
        k = KVector{1,Q,T}(ntuple(i -> i == n, Val(dimension(Q))))
        return :($var = $k)
    end
    printed = :(@info "Defined basis vectors: " *  join($vars, ", "))
    return Expr(:block, exs..., printed)
end

"""
    @basis_vars(Q, ::Type{T}; prefix = Metrics.blade_symbol(Q))

Generates variables in the global scope representing the 1-blade basis elements of `Q`.

!!! warning
    The default prefix for most algebras is 'e', which can cause problems with multiplying through
    juxtaposition: Julia interprets `2e0` as scientific notation. For the sake of clarity, use
    explicit multiplication (`2*e0`) or change the prefix.
    
    Lorentzian geometric algebras default to `γ` as the default prefix and do not have this problem.

# Examples
```julia-repl
julia> @basis_vars(VGA(3), prefix = :σ)

```
"""
macro basis_vars(Q, T, prefix)
    _Q = __module__.eval(Q)
    _T = __module__.eval(T)
    _prefix = __module__.eval(prefix)
    ex = basis_vars_expr(Val(_Q), _T, _prefix)
    return esc(:($ex))
end

macro basis_vars(Q, T)
    _Q = __module__.eval(Q)
    _T = __module__.eval(T)
    ex = basis_vars_expr(Val(_Q), _T)
    return esc(:($ex))
end

macro basis_vars(Q)
    _Q = __module__.eval(Q)
    ex = basis_vars_expr(Val(_Q))
    return esc(:($ex))
end

