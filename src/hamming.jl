"""
    CliffordNumbers.isevil(i::Integer) -> Bool

Determines whether a number is evil, meaning that its Hamming weight (sum of its binary digits) is
even.
"""
isevil(i::Integer) = iseven(count_ones(i))

"""
    CliffordNumbers.isodious(i::Integer) -> Bool

Determines whether a number is odious, meaning that its Hamming weight (sum of its binary digits) is
odd.
"""
isodious(i::Integer) = !isevil(i)

"""
    CliffordNumbers.evil_number(n::Integer)

Returns the nth evil number, with the first evil number (`n == 1`) defined to be 0.

Evil numbers are numbers which have an even Hamming weight (sum of its binary digits).
"""
evil_number(n::Integer) = (x = 2*(n-1); return x + !isevil(x))

"""
    CliffordNumbers.odious_number(n::Integer)

Returns the nth odious number, with the first odious number (`n == 1`) defined to be 1.

Odious numbers are numbers which have an odd Hamming weight (sum of its binary digits).
"""
odious_number(n::Integer) = (x = 2*(n-1); return x + !isodious(x))

"""
    CliffordNumbers.next_of_hamming_weight(n::Integer)

Returns the next integer with the same Hamming weight as `n`.
"""
function next_of_hamming_weight(n::Integer)
    c = n & -n
    r = n + c
    return div(xor(r, n) >>> 2, c) | r
end

# TODO: can this be done in a time under O(n) with constant space?
"""
    CliffordNumbers.hamming_number(w::Integer, n::Integer)

Gets the `n`th number with Hamming weight `w`. The first number with this Hamming weight (`n = 1`)
is `2^w - 1`.

# Example

```julia-repl
julia> CliffordNumbers.hamming_number(3, 2)
11
```
"""
function hamming_number(w::Integer, n::Integer)
    isone(w) && return 2^(n-1)
    result = 2^w - 1
    for _ in 1:(n-1)
        result = next_of_hamming_weight(result)
    end
    return result
end
