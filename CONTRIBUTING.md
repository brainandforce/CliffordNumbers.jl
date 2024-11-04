# Contributing to [CliffordNumbers.jl]

Thanks for taking an interest in contributing to this package!

## Prerequisites

I, the package author, believe that anyone who is willing to contribute to this package can do so.
However, you will need an understanding of both the Julia programming language and geometric algebra
in order to understand what this package is doing and where it can be improved.

For learning about Julia, the [official documentation][julia-docs] is a great starting point,
especially if you already know the fundamentals of programming.

For learning about geometric algebra, many resources exist. The implementation here was based off
of the one described in the textbook [Geometric Algebra for Computer Science][GA4CS], so I encourage
you to read this. If you prefer a more visual introduction to geometric algebra, there are many
excellent videos: most notably, [A Swift Introduction to Geometric Algebra][swift-introduction] by
sudgylacmoe.

And finally, you will need to know the basics of Git (and by extension, GitHub) to understand how to
contribute, especially if you wish to file a pull request.

If you're still not sure about how to contribute to this package, feel free to send me an email or
contact me on a different platform. While I can't guarantee my availabilty at convenient times for 
you or my capacity to help for extended periods, I will try my best to accomodate you: everyone
interested in programming deserves a warm welcome.

## Code style and formatting

**This code is written with the expectation that it will not just be used to perform calculations
with geometric algebras, but that it will be studied by people interested in learning about how 
geometric algebra is implemente in software.** Clarity of code is essential, and comments should be
used wherever something is unclear. If you had to take a significant amount of time to understand
what a piece of code is doing, please include explanatory comments, update the docstring, or rename
variables to clarify the code. I welcome PRs that serve only to clarify code.

At the moment, there isn't rigorous checking of the code format using external tools. Here are some 
rules to follow for the time being:
  * All indents are four spaces, like in most other Julia codebases.
  * Limit lines to 100 characters in all source files.
  * All modules, types, and data structures, as well as their aliases, should use `UpperCamelCase`.
  * All functions should use `lower_snake_case`.
  * All constants should use `ALLCAPS`, unless the constant is a type alias, in which case it should
    use `UpperCamelCase`.
  * All operators (such as `∧`) should be defined with their operator symbol and use aliases for
    verbose names (in this case, `wedge`). Any operators defined in this package must have a verbose
    name that uses ASCII characters.

## Interoperability

### Julia version

CliffordNumbers.jl will always support the latest Julia LTS at the time of release. At the time of 
writing, the LTS version is Julia 1.10.6. Support for earlier versions is provided if possible, but
it is not guaranteed.

Dropping support for any Julia version is considered a breaking change. At this stage of 
development, breaking changes are expected, and any changes that relies on language features that
are only available on the latest LTS are welcome. However, code using language features that would 
require the use of a version beyond the current LTS should not be included in the core codebase.
However, such code *is* acceptable in extensions for packages that only support a release after the
current LTS.

### Dependencies

At the time of writing, this package does not have any dependencies. Zero dependencies is not a
goal of the package: part of the reason for this is that I may spin some of the functionality of
this package off into separate subpackages (such as the `CliffordNumbers.Hamming` module for working
with Hamming weights of integers or other binary data).

However, this package is intended to be a relatively low-level component of the Julia ecosystem.
No dependencies on packages providing high-level functionality (such as GUIs) should ever be
included. In addition, this package has not needed to depend on packages that would ordinarily be
considered for a library like this, such as [LinearAlgebra], [StaticArrays], or [StaticArraysCore].
Inclusion of these libraries as dependencies will only be done if a pressing need for them can be
demonstrated.

### Extensions

If you foresee the need for interoperability with a different package that cannot be achieved by
simply loading the packages together, you are strongly encouraged to write a package extension.

As of the time of writing, this package is still in the 0.1 series, and supports Julia 1.8, which
does not have native support for package extensions. As of 0.2, support for Julia 1.8 will be
dropped, meaning that the native package extension functionality of the Julia package manager will
be utilized exclusively (this requires Julia 1.9).

#### Who should own the extension?

This package is MIT licensed, and if you consent to your contributions being released under the same
license or a fully compatible one, there is no problem with the package being hosted here. If you 
intend for the extension to be licensed under the same license as the other target package, it is
much simpler to host the extension with the target package.

If you wish to contribute an extension to this package under a free software license that is not
compatible with the MIT License, contact me so we can discuss how to approach this best.

Outside of that, here are some guidelines for where the package should be hosted:
  * I encourage you to host the extension with the other target package if you anticipate that the
    other target will work with or wrap types provided by CliffordNumbers.jl. This is especially
    true if the main issue you are trying to overcome is the fact that `AbstractCliffordNumber` does
    not subtype `AbstractArray`.
  * If your package extension intends to provide optimizations for performing operations in this
    package on different hardware (such as GPUs), the extension is likely better off hosted with
    CliffordNumbers.jl.

## Implementations of operations

Implementing operations in this package must prioritize correctness before anything else. While
performance considerations are critical, improvements to performance cannot come at the cost of
correctness or numerical stability.

### Testing and correctness

Testing of all functionality provided by this package is *absolutely essential* to ensure that all
results produced by this package are thoroughly correct. Like any other numerical library, any bug
that produces an incorrect answer for a given expression is a major issue that must be resolved
immediately for the sake of all dependents.

Any examples of expressions that generate incorrect results *must* be added to the test suite to
ensure that the same bug does not reappear in the code in a later version.

Although codecov can tell you which lines are missing test coverage, that does not mean testing of
that code is sufficient, and addition of more test cases is always welcome.

### Performance

This package is intended to be the fastest available implementation of arbitrary-dimensional
geometric algebras in the Julia ecosystem. To accomplish this, it heavily leverages tools provided
by Julia, such as generated functions, and functionality provided by hardware, such as SIMD 
registers and operations.

Any performance improvements in the implementations of operations are welcome. However, they must
work correctly on all platforms supported by Julia and they should not sacrifice correctness or
numerical stability for the sake of speed (unless the operation is implemented as a separate 
function, and documentation for the implementation *explicitly* states such a tradeoff is being 
made).

Performance regressions in a release are considered bugs and you are encouraged to file an issue
if you notice one. The only time a performance regression will not be fixed is if the drop in
performance is caused by the resolution of a correctness bug.

### Generated functions

CliffordNumbers.jl leverages [`@generated` functions][generated-functions] heavily to manually
unroll loops of fixed size and maximize performance for common operations, such as the various
products. However, they come at a cost: they are significantly more difficult to debug and must be
compiled for every possible combination of arguments.

In general, it is a good idea to use a generated function as a kernel for working with restricted
combinations of input types, and call that kernel with other methods. This is how products are
implemented: the generated methods of function `CliffordNumbers.mul` are only produced for cases
where the scalar types of the arguments are the same. A non-generated implementation which accepts
multivectors with different scalar types performs the promotion and then calls the generated
implementation.

[Optionally-generated functions][optionally-generated] (containing an `if @generated` statement) 
overcome some of the issues generated functions have. They are welcome as part of an implementation,
but only if their inclusion does not lead to significant performance regressions. Whether Julia
decides to use a generated or non-generated implementation is an implementation detail.

### SIMD

Understanding how SIMD operations work can help you understand why code in this package is written
the way it is, and how you may be able to improve the performance of slower operations in this
package. We encourage you to learn how code lowering tools like `@code_typed`, `@code_llvm` and
`@code_native` work if you want to work on performance optimizations.

It is okay to rely on hardware details, such as SIMD register width, to generate optimized code.
However, this information must be acquired in an OS-independent manner, and reasonable fallbacks 
must exist in case that information cannot be obtained.

## Documentation

Thorough documentation is critical for ordinary users of this package and for future developers. As
mentioned before, the goals of CliffordNumbers.jl include writing code that is pedagogically useful
as well as correct and performant.

**All symbols must have a docstring, regardless of whether or not they are exported, public, or 
private.** The format of the docstring should follow that of Julia Base: the name of the function 
must be given on the first line, indented by four spaces to trigger monospace code formatting.

Not all methods of a function need to be documented if their semantics are identical to other
methods. However, function methods with significantly different semantics should have separate
docstrings.

Unexported symbols should be prefixed with the module name. In Julia 1.11 or later, unexported
names are treated as private, and their docstrings will include a warning indicating that the symbol
is a private method, so it is not necessary to include this information in the docstring.

### Documenter

This package includes online documentation hosted by GitHub and powered by [Documenter]. Building of
the documentation is part of the continuous integration pipeline, and pull requests will only be
merged when the documenation builds without error.

One of the most common causes of documentation build failures is that all docstrings are expected to
be included somewhere in the documentation exactly once. If you documented a symbol -- a new type,
function, or even a method of an already existing function -- you must include the docstring in an 
appropriate API reference page.

Including a docstring more than once in all of the documentation will also trigger an error. If you 
wish to include a docstring more than once, you must include the parameter `canonical=false` in the 
code block. The canonical docstring must be included in the API reference.

Outside of API documentation and docstring inclusion, Documenter supports the use of LaTeX, and you
are strongly encouraged to rely on it wherever mathematical notation is needed.

## Licensing

This software falls under the [MIT License](LICENSE), and all contributions to the core codebase
(not including extensions) must be released under the MIT License or a compatible license.

For extensions, licensing issues may be resolved by contributing the extension to the other target
package, if the intended license is identical to that of the other package. However, I am not 
opposed to hosting extensions that are under different licenses, though this will require further
discussion.

I will *never* merge any pull requests that include code that falls under a nonfree software license
per the [GNU Free Software Defintion][GNU FSD], [Debian Free Software Guidelines][DFSG], or the 
[Open Source Definition][OSD]. They will be closed immediately.

[CliffordNumbers.jl]:   https://github.com/brainandforce/CliffordNumbers.jl
[julia-docs]:           https://docs.julialang.org/
[GA4CS]:                https://geometricalgebra.org/
[swift-introduction]:   https://www.youtube.com/watch?v=60z_hpEAtD8
[generated-functions]:  https://docs.julialang.org/en/v1/manual/metaprogramming/#Generated-functions
[optionally-generated]: https://docs.julialang.org/en/v1/manual/metaprogramming/#Optionally-generated-functions
[LinearAlgebra]:        https://github.com/JuliaLang/LinearAlgebra.jl
[StaticArrays]:         https://github.com/JuliaArrays/StaticArrays.jl
[StaticArraysCore]:     https://github.com/JuliaArrays/StaticArraysCore.jl
[Documenter]:           https://github.com/JuliaDocs/Documenter.jl
[GNU FSD]:              https://www.gnu.org/philosophy/free-sw.en.html
[DFSG]:                 https://wiki.debian.org/DebianFreeSoftwareGuidelines
[OSD]:                  https://opensource.org/osd
