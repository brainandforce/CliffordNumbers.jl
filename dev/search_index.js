var documenterSearchIndex = {"docs":
[{"location":"getting_started/#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"getting_started/#The-CliffordNumber-data-type","page":"Getting started","title":"The CliffordNumber data type","text":"","category":"section"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"A CliffordNumber{Q,T,L} is a Clifford number associated with a QuadraticForm Q, a backing Real or Complex type T, and a length L. The length parameter is redundant, and in many cases, it may be omitted without consequence.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"!!! warning Although the length may be omitted in many cases, it's important to remember that a CliffordNumber{Q,T} is not a concrete type. This is important when creating an Array or other container of CliffordNumber elements.","category":"page"},{"location":"getting_started/#Internals","page":"Getting started","title":"Internals","text":"","category":"section"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"A CliffordNumber{Q,T,L} is backed by an NTuple{L,T} where T<:Union{Real,Complex}. The coefficients, however, are not indexed in grade order as is done canonically in most resources.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"!!! danger Read that again: CliffordNumber indexing is not done in grade order.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"Instead, the coefficients are arranged in a binary counted fashion, which allows for better SIMD optimization.","category":"page"},{"location":"getting_started/#Constructing-a-Clifford-number","page":"Getting started","title":"Constructing a Clifford number","text":"","category":"section"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"The inner constructor for CliffordNumber is CliffordNumber{Cl,T,L}(x), where x is any type that can be converted to an NTuple{L,T}. However, in many cases, the type parameters are redundant, particularly L. For this reason, more constructors exist.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"In general, one can use a Vararg constructor to directly input the values.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"julia> CliffordNumber{APS}(1, 2, 3, 4, 5, 6, 7, 8)\nCliffordNumber{APS,Int}(1, 2, 3, 4, 5, 6, 7, 8)","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"Clifford numbers may also be constructed from real numbers, generating a scalar-valued CliffordNumber:","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"julia> CliffordNumber{APS}(1)\nCliffordNumber{APS,Int}(1, 0, 0, 0, 0, 0, 0, 0)","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"When constructing a CliffordNumber from complex numbers, the type parameters become more important. By default, it is assumed that the element type of a CliffordNumber is a Real. If a complex CliffordNumber is desired, this must be stated explicitly.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"julia> CliffordNumber{APS}(1 + im)\nCliffordNumber{APS,Int}(1, 0, 0, 0, 0, 0, 0, 1)\n\njulia> CliffordNumber{APS,Complex}(1 + im)\nCliffordNumber{APS,Complex{Int}}(1 + im, 0, 0, 0, 0, 0, 0, 0)","category":"page"},{"location":"getting_started/#Quadratic-forms","page":"Getting started","title":"Quadratic forms","text":"","category":"section"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"Before getting started with Clifford numbers, it's important to understand how the dimensionality of the space is stored. Unlike with other data types such as StaticArrays.jl's SVector, the total number of dimensions in the space is not all the information that needs to be stored. Each basis vector of the space may square to a positive number, negative number, or zero, defining the quadratic form associated with the Clifford algebra. This information needs to be tracked as a type parameter for CliffordNumber.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"To handle this, the QuadraticForm{P,Q,R} type is used to store information about the quadratic form. In this type, P represents","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"!!! note By convention, the QuadraticForm type is not instantiated when used as a type parameter for CliffordNumber instances.","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"CliffordNumbers.jl provides the following aliases for common algebras:","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"| Algebra    | Alias                  | Note                                           | | VGA{D}   | QuadraticForm{D,0,0} | Vanilla/vector geometric algebra               | | PGA{D}   | QuadraticForm{D,0,1} | Projective geometric algebra                   | | APS      | QuadraticForm{3,0,0} | Algebra of physical space                      | | STA      | QuadraticForm{1,3,0} | Spacetime algebra. By default, uses a -+++     | |            |                        | convention to distinguish it from a conformal  | |            |                        | geometric algebra.                             |","category":"page"},{"location":"getting_started/","page":"Getting started","title":"Getting started","text":"Currently, an alias for conformal geometric algebras (CGA{D}) does not exist, as it requires some type parameter trickery that hasn't been figured out yet.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = CliffordNumbers","category":"page"},{"location":"#CliffordNumbers","page":"Home","title":"CliffordNumbers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CliffordNumbers.jl is a package that provides fully static multivectors (Clifford numbers) in arbitrary dimensions and metrics. While in many cases, sparse representations of multivectors are more efficient, for spaces of low dimension, dense static representations may provide a performance and convenience advantage.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CliffordNumbers]","category":"page"},{"location":"#CliffordNumbers.BaseNumber","page":"Home","title":"CliffordNumbers.BaseNumber","text":"CliffordNumbers.BaseNumber\n\nUnion of subtypes of Number provided in the Julia Base module: Real and Complex. This encompasses all types that may be used to construct a CliffordNumber.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.APS","page":"Home","title":"CliffordNumbers.APS","text":"APS\n\nThe algebra of physical space, Cl(3,0,0). An alias for QuadraticForm{3,0,0}.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.AbstractCliffordNumber","page":"Home","title":"CliffordNumbers.AbstractCliffordNumber","text":"AbstractCliffordNumber{Q,T}\n\nAn element of a Clifford algebra, often referred to as a multivector, with quadratic form Q and element type T.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.BitIndex","page":"Home","title":"CliffordNumbers.BitIndex","text":"BitIndex{Q<:QuadraticForm}\n\nA representation of an index corresponding to a basis blade of the geometric algebra with quadratic form Q.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.BitIndex-Union{Tuple{Vararg{Integer}}, Tuple{Q}} where Q","page":"Home","title":"CliffordNumbers.BitIndex","text":"BitIndex{Q}(x::Integer...)\nBitIndex{Q}(v::AbstractVector{<:Integer})\n\nConstructs a BitIndex{Q} from a list of integers that represent the basis vectors of the space.\n\nThis package uses a lexicographic convention for basis blades: in the algebra of physical space, the basis bivectors are {e₁e₂, e₁e₃, e₂e₃}. The sign of the BitIndex{Q} is negative when the parity of the basis vector permutation is odd.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.BitIndices","page":"Home","title":"CliffordNumbers.BitIndices","text":"BitIndices{Q<:QuadraticForm}\n\nRepresents a range of valid BitIndex objects for a given quadratic form.\n\nFor sparse representations, such as KVector{K...}, this is not the most efficient way to iterate through all elements, as it includes indices that are known to be zero.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.BitIndices-Union{Tuple{Type{<:AbstractCliffordNumber{Q}}}, Tuple{Q}} where Q","page":"Home","title":"CliffordNumbers.BitIndices","text":"BitIndices(x::AbstractCliffordNumber{Q}) -> BitIndices{Q}()\nBitIndices(T::Type{<:AbstractCliffordNumber{Q}}) -> BitIndices{Q}()\n\nConstructs a BitIndices object associated with a Clifford number or its type.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.CliffordNumber","page":"Home","title":"CliffordNumbers.CliffordNumber","text":"CliffordNumber{Q,T,L}\n\nA dense multivector (or Clifford number), with quadratic form Q, element type T, and length L (which depends entirely on Q).\n\nThe coefficients are ordered by taking advantage of the natural binary structure of the basis. The grade of an element is given by the Hamming weight of its index. For the algebra of physical space, the order is: 1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃ = i. This order allows for more aggressive SIMD optimization when calculating the geometric product.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.KVector","page":"Home","title":"CliffordNumbers.KVector","text":"KVector{K,Q,T,L}\n\nA multivector consisting only linear combinations of basis blades of grade K - in other words, a k-vector.\n\nk-vectors have binomial(dimension(Q), K) components.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.QuadraticForm","page":"Home","title":"CliffordNumbers.QuadraticForm","text":"CliffordNumbers.QuadratricForm\n\nRepresents a quadratric with P dimensions which square to +1, Q dimensions which square to -1, and R dimensions which square to 0, in that order.\n\nBy convention, this type is used as a tag, and is never instantiated.\n\n\n\n\n\n","category":"type"},{"location":"#CliffordNumbers.STA","page":"Home","title":"CliffordNumbers.STA","text":"STA\n\nSpacetime algebra with a mostly negative signature (particle physicist's convention), Cl(1,3,0). An alias for QuadraticForm{1,3,0}.\n\nThe negative signature is used by default to distinguish this algebra from conformal geometric algebras, which use a mostly positive signature by convention.\n\n\n\n\n\n","category":"type"},{"location":"#Base.:*-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"Base.:*","text":"*(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the geometric product between multivectors/Clifford numbers x and y which share the quadratic form Q.\n\n\n\n\n\n","category":"method"},{"location":"#Base.:~-Tuple{CliffordNumber}","page":"Home","title":"Base.:~","text":"reverse(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n~(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n\nCalculate the reverse of a Clifford number. This effectively reverses the products that form the basis blades, or in other words, reverses the order of the geometric product that resulted in x.\n\n\n\n\n\n","category":"method"},{"location":"#Base.abs-Tuple{CliffordNumber}","page":"Home","title":"Base.abs","text":"abs2(x::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the norm of x, equal to sqrt(scalar_product(x, ~x)).\n\n\n\n\n\n","category":"method"},{"location":"#Base.abs2-Tuple{CliffordNumber}","page":"Home","title":"Base.abs2","text":"abs2(x::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the squared norm of x, equal to scalar_product(x, ~x).\n\n\n\n\n\n","category":"method"},{"location":"#Base.conj-Tuple{CliffordNumber}","page":"Home","title":"Base.conj","text":"conj(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n\nCalculates the Clifford conjugate of a Clifford number x. This is equal to grade_involution(reverse(x)).\n\n\n\n\n\n","category":"method"},{"location":"#Base.exp-Tuple{CliffordNumber}","page":"Home","title":"Base.exp","text":"exp(x::CliffordNumber{Q}) -> CliffordNumber{Q,<:AbstractFloat}\n\nReturns the natural exponential of a Clifford number.\n\nFor special cases where m squares to a scalar, the following shortcuts can be used to calculate exp(x):\n\nWhen x^2 < 0: exp(x) === cos(abs(x)) + x * sin(abs(x)) / abs(x)\nWhen x^2 > 0: exp(x) === cosh(abs(x)) + x * sinh(abs(x)) / abs(x)\nWhen x^2 === 0: exp(x) == 1 + x\n\nSee also: exppi, exptau.\n\n\n\n\n\n","category":"method"},{"location":"#Base.real-Tuple{CliffordNumber{<:QuadraticForm, <:Real}}","page":"Home","title":"Base.real","text":"real(x::CliffordNumber{Q,T<:Real}) = T\n\nReturn the real (scalar) portion of a real Clifford number. \n\n\n\n\n\n","category":"method"},{"location":"#Base.reverse-Tuple{CliffordNumber}","page":"Home","title":"Base.reverse","text":"reverse(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n~(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n\nCalculate the reverse of Clifford number x. This effectively reverses the products that form the basis blades, or in other words, reverses the order of the geometric product that resulted in x.\n\n\n\n\n\n","category":"method"},{"location":"#Base.sign-Union{Tuple{R}, Tuple{Q}, Tuple{P}, Tuple{Type{QuadraticForm{P, Q, R}}, Integer}} where {P, Q, R}","page":"Home","title":"Base.sign","text":"sign(::Type{QuadraticForm{P,Q,R}}, i::Integer) -> Int8\n\nGets the sign associated with dimension i of a quadratric form.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.CGA-Tuple{Any}","page":"Home","title":"CliffordNumbers.CGA","text":"CGA(D) -> Type{QuadraticForm{D+1,1,0}}\n\nCreates the type of a quadratic form associated with a conformal geometric algebra (CGA) of dimension D.\n\nFor reasons of type stability, avoid calling this function without constant arguments.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.PGA-Tuple{Any}","page":"Home","title":"CliffordNumbers.PGA","text":"PGA(D) -> Type{QuadraticForm{D,0,1}}\n\nCreates the type of a quadratic form associated with a projective geometric algebra (PGA) of dimension D.\n\nFor reasons of type stability, avoid calling this function without constant arguments.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.VGA-Tuple{Any}","page":"Home","title":"CliffordNumbers.VGA","text":"VGA(D) -> Type{QuadraticForm{D,0,0}}\n\nCreates the type of a quadratic form associated with a vector/vanilla geometric algebra (VGA) of dimension D.\n\nFor reasons of type stability, avoid calling this function without constant arguments.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers._sort_with_parity!-Tuple{AbstractVector{<:Real}}","page":"Home","title":"CliffordNumbers._sort_with_parity!","text":"CliffordNumbers._sort_with_parity!(v::AbstractVector{<:Real}) -> Tuple{typeof(v),Bool}\n\nPerforms a parity-tracking insertion sort of v, which modifies v in place. The function returns a tuple containing v and the parity, which is true for an odd permutation and false for an even permutation. This is implemented with a modified insertion sort algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.dot-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.dot","text":"dot(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the dot product of x and y.\n\nFor basis blades A of grade m and B of grade n, the dot product is equal to the left contraction when m >= n and is equal to the right contraction when n >= m.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.dual-Union{Tuple{CliffordNumber{Q}}, Tuple{Q}} where Q","page":"Home","title":"CliffordNumbers.dual","text":"dual(x::CliffordNumber) -> CliffordNumber\n\nCalculates the dual of x, which is equal to the left contraction of x with the inverse of the pseudoscalar. However, \n\nNote that the dual has some properties that depend on the dimension and quadratic form:\n\nThe inverse of the unit pseudoscalar depends on the dimension of the space. Therefore, the\n\nperiodicity of \n\nIf the metric is degenerate, the dual is not unique.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.elementwise_product-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}, BitIndex{Q}, BitIndex{Q}}} where Q","page":"Home","title":"CliffordNumbers.elementwise_product","text":"CliffordAlgebra.elementwise_product(\n    x::CliffordNumber{Q},\n    y::CliffordNumber{Q},\n    a::BitIndex{Q},\n    b::BitIndex{Q}\n)\n\nCalculates the geometric product between the element of x indexed by a and the element of y indexed by b.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.evil_number-Tuple{Integer}","page":"Home","title":"CliffordNumbers.evil_number","text":"CliffordNumbers.evil_number(n::Integer)\n\nReturns the nth evil number, with the first evil number (n == 1) defined to be 0.\n\nEvil numbers are numbers which have an even Hamming weight (sum of its binary digits).\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.exppi-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.exppi","text":"exppi(x::CliffordNumber)\n\nReturns the natural exponential of π * x with greater accuracy than exp(π * x) in the case where x^2 is a negative scalar.\n\nSee also: exp, exptau.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.exptau-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.exptau","text":"exptau(x::CliffordNumber)\n\nReturns the natural exponential of 2π * x with greater accuracy than exp(2π * x) in the case where x^2 is a negative scalar.\n\nSee also: exp, exppi.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.grade_involution-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.grade_involution","text":"grade_involution(x::CliffordNumber{Q,T}) -> CliffordNumber{Q,T}\n\nCalculates the grade involution of Clifford number x. This effectively multiplies all of the basis vectors of the space by -1, which makes elements of odd grade flip sign.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.hestenes_product-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.hestenes_product","text":"hestenes_product(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nReturns the Hestenes product: this is equal to the dot product given by dot(x, y) but is equal to to zero when either x or y is a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.isevil-Tuple{Integer}","page":"Home","title":"CliffordNumbers.isevil","text":"CliffordNumbers.isevil(i::Integer) -> Bool\n\nDetermines whether a number is evil, meaning that its Hamming weight (sum of its binary digits) is even.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.isodious-Tuple{Integer}","page":"Home","title":"CliffordNumbers.isodious","text":"CliffordNumbers.isodious(i::Integer) -> Bool\n\nDetermines whether a number is odious, meaning that its Hamming weight (sum of its binary digits) is odd.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.ispseudoscalar-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.ispseudoscalar","text":"ispseudoscalar(m::CliffordNumber)\n\nDetermines whether the Clifford number m is a pseudoscalar, meaning that it has no components with grades below the one equal to the dimension of the space.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.isscalar-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.isscalar","text":"isscalar(m::CliffordNumber)\n\nDetermines whether the Clifford number m is a scalar, meaning that it has no components with grades above zero.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.left_contraction-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.left_contraction","text":"left_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the left contraction of x and y.\n\nFor basis blades A of grade m and B of grade n, the left contraction is zero if n < m, otherwise it is grade_select(A*B, n-m).\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.normalize-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.normalize","text":"normalize(x::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nNormalizes x so that its magnitude (as calculated by abs2(x)) is 1.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.odious_number-Tuple{Integer}","page":"Home","title":"CliffordNumbers.odious_number","text":"CliffordNumbers.odious_number(n::Integer)\n\nReturns the nth odious number, with the first odious number (n == 1) defined to be 1.\n\nOdious numbers are numbers which have an odd Hamming weight (sum of its binary digits).\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.right_contraction-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.right_contraction","text":"right_contraction(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the right contraction of x and y.\n\nFor basis blades A of grade m and B of grade n, the right contraction is zero if m < n, otherwise it is grade_select(A*B, m-n).\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.sandwich-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.sandwich","text":"sandwich(x::CliffordNumber{Q}, y::CliffordNumber{Q})\n\nCalculates the sandwich product of x with y: ~y * x * y, but with corrections for numerical stability. \n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.scalar_product-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.scalar_product","text":"scalar_product(x::CliffordNumber{Q,T1}, y::CliffordNumber{Q,T2}) -> promote_type(T1,T2)\n\nCalculates the scalar product of two Clifford numbers with quadratic form Q. The result is a Real or Complex number. This can be converted back to a CliffordNumber.\n\nThis is equal to grade_select(x*y, 0) but is significantly more efficient.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.select_grade-Tuple{CliffordNumber, Integer}","page":"Home","title":"CliffordNumbers.select_grade","text":"select_grade(x::CliffordNumber, g::Integer)\n\nReturns a multivector similar to x where all elements not of grade g are equal to zero.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.sign_of_mult-Union{Tuple{R}, Tuple{Q}, Tuple{P}, Tuple{BitIndex{QuadraticForm{P, Q, R}}, BitIndex{QuadraticForm{P, Q, R}}}} where {P, Q, R}","page":"Home","title":"CliffordNumbers.sign_of_mult","text":"CliffordNumbers.sign_of_mult(a::T, b::T) where T<:BitIndex{QuadraticForm{P,Q,R}} -> Int8\n\nReturns an Int8 that carries the sign associated with the multiplication of two basis blades of Clifford/geometric algebras of the same quadratic form.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.signbit_of_mult-Tuple{Unsigned, Unsigned}","page":"Home","title":"CliffordNumbers.signbit_of_mult","text":"CliffordNumbers.signbit_of_mult(a::Integer, [b::Integer]) -> Bool\nCliffordNumbers.signbit_of_mult(a::BitIndex, [b::BitIndex]) -> Bool\n\nCalculates the sign bit associated with multiplying basis elements indexed with bit indices supplied as either integers or BitIndex instances. The sign bit flips when the order of a and b are reversed, unless a === b. \n\nAs with Base.signbit(), true represents a negative sign and false a positive sign. However, in degenerate metrics (such as those of projective geometric algebras) the sign bit may be irrelevant as the multiplication of those basis blades would result in zero.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.signmask","page":"Home","title":"CliffordNumbers.signmask","text":"CliffordNumbers.signmask(T::Type{<:Integer}, signbit::Bool = true) -> T1\n\nGenerates a signmask, or a string of bits where the only 1 bit is the sign bit. If signbit is set to false, this returns a string of bits.\n\n\n\n\n\n","category":"function"},{"location":"#CliffordNumbers.subscript_string-Tuple{Number}","page":"Home","title":"CliffordNumbers.subscript_string","text":"CliffordNumbers.subscript_string(x::Number) -> String\n\nProduces a string representation of a number in subscript format.\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.undual-Union{Tuple{CliffordNumber{Q}}, Tuple{Q}} where Q","page":"Home","title":"CliffordNumbers.undual","text":"undual(x::CliffordNumber) -> CliffordNumber\n\nCalculates the undual of x, which is equal to the left contraction of x with the pseudoscalar. This function can be used to reverse the behavior of dual().\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.versor_inverse-Tuple{CliffordNumber}","page":"Home","title":"CliffordNumbers.versor_inverse","text":"versor_inverse(x::CliffordNumber)\n\nCalculates the versor inverse of x, equal to x / scalar_product(x, ~x), so that x * inv(x) == inv(x) * x == 1.\n\nThe versor inverse is only guaranteed to be an inverse for blades and versors. Not all Clifford numbers have a well-defined inverse, since Clifford numbers have zero divisors (for instance, in the algebra of physical space, 1 + e₁ has a zero divisor).\n\n\n\n\n\n","category":"method"},{"location":"#CliffordNumbers.wedge-Union{Tuple{Q}, Tuple{CliffordNumber{Q}, CliffordNumber{Q}}} where Q","page":"Home","title":"CliffordNumbers.wedge","text":"wedge(x::CliffordNumber{Q}, y::CliffordNumber{Q}) -> CliffordNumber{Q}\n\nCalculates the wedge (outer) product of two Clifford numbers x and y with quadratic form Q.\n\n\n\n\n\n","category":"method"}]
}
