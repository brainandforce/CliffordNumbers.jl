# Entry point for the CliffordNumbers benchmark suite.
#
#   julia --project=bench bench/runbenchmarks.jl
#   BENCH_SECONDS=5 julia -t auto --project=bench bench/runbenchmarks.jl
#
# Activating and instantiating `bench/` is also done here so the script runs
# against the in-repo CliffordNumbers (pinned via the `[sources]` entry) with a
# single `julia bench/runbenchmarks.jl` invocation.

using Pkg
if abspath(Base.active_project()) != abspath(joinpath(@__DIR__, "Project.toml"))
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
end

include("harness.jl")
include("bench_complex.jl")
include("bench_quaternion.jl")
include("bench_multivector.jl")

function main()
    failures = 0
    failures += report("bench_complex: even VGA(2) ≅ ℂ", complex_benchmarks())
    failures += report("bench_quaternion: even VGA(3) ≅ ℍ", quaternion_benchmarks())
    failures += report("bench_kvector: geometric primitives & incidence (compact vs dense)", kvector_benchmarks())
    failures += report("bench_even: rotors & motors (EvenCliffordNumber vs dense)", even_benchmarks())
    failures += report("bench_odd: reflections (OddCliffordNumber vs dense)", odd_benchmarks())
    failures += report("bench_general: full multivectors & inverses (CliffordNumber)", general_benchmarks())
    failures += report("bench_simd: batched broadcast at increasing N (compact vs dense)", simd_benchmarks())
    println()
    if failures == 0
        println("All algebra-drift checks passed.")
    else
        println(failures, " workload(s) FAILED the algebra-drift check.")
    end
    return failures
end

failures = main()
exit(failures == 0 ? 0 : 1)
