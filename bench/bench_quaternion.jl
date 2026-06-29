# Rotor (≅ ℍ) workloads: the batched rows targeted by PR-ArrayVectorization.
#
# `Rotor{T} = EvenCliffordNumber{VGA(3),T,4}` is exactly the quaternion shape.
# Compose / chain-compose are checked against ℍ multiplication through the
# `as_hamilton` ring homomorphism; slerp is checked through CliffordNumbers' own
# `Quaternion` conversion (matching the package test); the sandwich is checked by
# norm preservation (and a vanishing trivector part), which is convention-free.

# A grade-1 result preserves the input norm and carries no trivector part.
sandwich_ok(out, vnorm2) =
    (g = Tuple(out); isapprox(g[1]^2 + g[2]^2 + g[3]^2, vnorm2; rtol = 1e-10) && abs(g[4]) < 1e-9)

function quaternion_benchmarks()
    rng = MersenneTwister(0x0707_50A7)
    T = Float64
    M = 1000           # batch size for the per-op (single-call) rows
    results = WResult[]

    # --- compose, per op ---
    rs1 = [randrotor(rng) for _ in 1:M]; rs2 = [randrotor(rng) for _ in 1:M]
    qs1 = as_hamilton.(rs1); qs2 = as_hamilton.(rs2)
    push!(results, run_per_op("compose (per op)", M,
        () -> map(*, rs1, rs2),
        () -> map(*, qs1, qs2),
        (c, ref) -> approxall(as_hamilton.(c), ref; rtol = sqrt(eps(T)))))

    # --- single-vector sandwich vs quaternion rotation of a 3-tuple, per op ---
    rsv = randrotor(rng); qsv = Quaternion(rsv)
    vsv = [Vec3{T}(ntuple(_ -> randn(rng), 3)) for _ in 1:M]
    vsvt = [Tuple(x)[1:3] for x in vsv]
    n2sv = abs2.(vsv)
    push!(results, run_per_op("sandwich (per op)", M,
        () -> map(v -> ga_rotate(rsv, v), vsv),
        () -> map(vt -> q_rotate(qsv, vt), vsvt),
        (outs, _) -> all(((o, n2),) -> sandwich_ok(o, n2), zip(outs, n2sv))))

    # --- batched sandwich at increasing N ---
    for (N, label) in ((100, "N=100"), (100_000, "N=100k"), (1_000_000, "N=1M"))
        vs = [Vec3{T}(ntuple(_ -> randn(rng), 3)) for _ in 1:N]
        vts = [Tuple(x)[1:3] for x in vs]
        n2 = abs2.(vs)
        rb = randrotor(rng); qb = Quaternion(rb)
        push!(results, run_workload("batched sandwich ($label)",
            () -> ga_rotate.(Ref(rb), vs),
            () -> q_rotate.(Ref(qb), vts),
            (outs, _) -> all(((o, nn),) -> sandwich_ok(o, nn), zip(outs, n2))))
    end

    # --- slerp, per op (checked through CliffordNumbers' own conversion) ---
    as = [randrotor(rng) for _ in 1:M]; bs = [randrotor(rng) for _ in 1:M]
    qas = Quaternion.(as); qbs = Quaternion.(bs)
    push!(results, run_per_op("slerp (per op)", M,
        () -> map((a, b) -> slerp(a, b, 0.5), as, bs),
        () -> map((qa, qb) -> slerp(qa, qb, 0.5), qas, qbs),
        (c, ref) -> approxall(Quaternion.(c), ref; rtol = sqrt(eps(T)))))

    # --- chain-compose, N rotors reduced left-to-right ---
    NC = 10_000
    rs = [randrotor(rng) for _ in 1:NC]
    qs = as_hamilton.(rs)
    push!(results, run_workload("chain-compose (N=10k)",
        () -> reduce(*, rs),
        () -> reduce(*, qs),
        (c, ref) -> isapprox(as_hamilton(c), ref; rtol = 1e-6, atol = 1e-8)))

    return results
end
