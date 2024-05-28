#! /usr/bin/env -S julia -O3

using CliffordNumbers
using StaticArrays
using BenchmarkTools

function run_benchmarks()
    # generate random APS multivectors
    c1 = CliffordNumber{VGA(3)}(ntuple(_ -> rand(Float32), Val(8)))
    c2 = CliffordNumber{VGA(3)}(ntuple(_ -> rand(Float32), Val(8)))

    # generate random 2×2 complex matrices
    m1 = SMatrix{2,2}(ntuple(_ -> rand(Complex{Float32}), Val(4)))
    m2 = SMatrix{2,2}(ntuple(_ -> rand(Complex{Float32}), Val(4)))

    # Tuned benchmarks
    bm = tune!(@benchmarkable($m1 * $m2))
    bc = tune!(@benchmarkable($c1 * $c2))

    return (run(bm), run(bc))
end

function julia_main()::Cint
    (bm, bc) = run_benchmarks()
    println(stdout, "Results for SMatrix{2,2,Complex{Float32}} multiply:")
    display(bm)
    println(stdout, "Results for CliffordNumber{VGA(3),Float32} multiply:")
    display(bc)
    return 0
end

julia_main()
