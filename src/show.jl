#---Pretty printing methods for exported types-----------------------------------------------------#
"""
    CliffordNumbers.subscript_string(x::Number) -> String

Produces a string representation of a number in subscript format.
"""
function subscript_string(x::Number)
    str = collect(string(x))
    for (n,c) in enumerate(str)
        ('0' <= c <= '9') && (str[n] = c + 0x2050)
        (c === '-') && (str[n] = '₋')
        (c === '+') && (str[n] = '₊')
    end
    return String(str)
end

"""
    CliffordNumbers.to_basis_str(b::BitIndex; [label], [pseudoscalar])

Creates a string representation of the basis element given by `b`.

The `label` parameter determines the symbol used for every basis element. In general, this defaults
to `e`, but a few special cases use different symbols:
  * For the algebra of physical space, `σ` is used.
  * For the spacetime algebras of either signature, `γ` is used.

If `pseudoscalar` is set, the pseudoscalar may be printed using a different symbol from the rest of
the basis elements.
"""
function to_basis_str(b::BitIndex{Q}; label = nothing, pseudoscalar = nothing) where Q
    if iszero(grade(b))
        return ""
    elseif grade(b) === dimension(Q) && !(isnothing(pseudoscalar) || isempty(pseudoscalar))
        return string(pseudoscalar)
    end
    if isnothing(label) || isempty(label)
        if Q === QuadraticForm{3,0,0}
            label = "σ"
        elseif Q === QuadraticForm{3,1,0} || QuadraticForm{1,3,0}
            label = "γ"
        else
            label = "e"
        end
    end
    return join(
        (label * subscript_string(n+1))^!iszero(2^n & UInt(b)) for n in 0:(dimension(Q) - 1)
    )
end

# Generic pretty-print method for all AbstractCliffordNumber instances
function show(io::IO, ::MIME"text/plain", x::AbstractCliffordNumber{Q}) where Q
    summary(io, x)
    # Flag to mark when we've *found the first nonzero* element
    ffn = false
    # Print the scalar component first
    if !iszero(x[BitIndex{Q}()])
        print(io, x[BitIndex{Q}()])
        ffn = true
    end
    # Loop through all the grades
    for n in 1:dimension(Q)
        # Group together all indices corresponding to the given grade
        for i in Iterators.filter(i -> grade(i) == n, BitIndices(x))
            if !iszero(x[i])
                print(
                    io, " "^ffn, (sign(x[i]) > 0 ? "+"^ffn : "-"), " "^ffn, abs(x[i]), 
                    "*"^(x[i] isa Bool), to_basis_str(i)
                )
                ffn = true
            end
        end
    end
    # If we got through the whole process, just print the zero
    # We reference the actual element so we can print signed zero if needed
    !ffn && print(io, x[BitIndex{Q}()])
end
