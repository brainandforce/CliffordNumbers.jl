# Multivector workloads: KVector / EvenCliffordNumber / OddCliffordNumber / CliffordNumber.
#
# Unlike the ℂ/ℍ benches, the general Clifford types have no single stdlib primitive
# to compare against. Instead each workload times the compact representation a real
# application reaches for (`KVector` for geometric primitives, `EvenCliffordNumber`
# for rotors/motors, `OddCliffordNumber` for reflections) against the identical
# computation carried out on the dense `CliffordNumber`. The ratio is the payoff of
# choosing the compact type, and the drift check is convention-free coefficient
# equality between the two representations (`dense`/`densematch` live in
# `harness.jl`).
#
# The grades, algebras, and exact operator pipelines mirror those used by
# ComputationalGeometricAlgebra.jl (PGA/CGA/VGA `meet`/`join`/`project`/`move`/
# `reflect`/`rotate`): e.g. in PGA(3) a point is a grade-3 `KVector`, a line is
# grade-2, a plane is grade-1, and a motor is an `EvenCliffordNumber`; `meet` is the
# wedge `∧`, `join` is the regressive `∨`, and a rigid motion is the sandwich
# `M X M̃`.

const RT = sqrt(eps(Float64))   # drift-check tolerance for the compact-vs-dense comparison

# A random grade-`K` k-vector of algebra `Q` with `Float64` scalars.
randk(rng::AbstractRNG, K::Integer, Q) = KVector{K,Q,Float64}(ntuple(_ -> randn(rng), binomial(dimension(Q), K)))

# A random even versor (rotor / motor) of algebra `Q`, built as `exp(B)` of a random
# bivector, as Rotation/Motor representations are generated in applications.
randversor(rng::AbstractRNG, Q) = exp(randk(rng, 2, Q))

# A random even versor as the geometric product of two 1-vectors. Used for CGA, whose
# mixed signature can route `exp` of a random bivector through the Taylor path; the
# two-reflection product is a genuine conformal versor and stays on the closed path.
randeven2(rng::AbstractRNG, Q) = randk(rng, 1, Q) * randk(rng, 1, Q)

# A random odd versor of VGA(3): the geometric product of three 1-vectors (an
# improper isometry, e.g. a reflection composed with a rotation).
randodd(rng::AbstractRNG, Q) = randk(rng, 1, Q) * randk(rng, 1, Q) * randk(rng, 1, Q)

# The sandwich/versor product `R X R̃`: the rotation / rigid-motion / reflection kernel.
sandwich(R, x) = R * x * reverse(R)

#---Per-op registration helpers (compact vs dense, reported per element)----------------------------#

function cvd1!(rs, name, M, op, xs)
    dxs = dense.(xs)
    push!(rs, run_per_op(name, M, () -> map(op, xs), () -> map(op, dxs),
        (c, r) -> densematch(c, r; rtol = RT)))
end

function cvd2!(rs, name, M, op, xs, ys)
    dxs = dense.(xs); dys = dense.(ys)
    push!(rs, run_per_op(name, M, () -> map(op, xs, ys), () -> map(op, dxs, dys),
        (c, r) -> densematch(c, r; rtol = RT)))
end

function cvd3!(rs, name, M, op, xs, ys, zs)
    dxs = dense.(xs); dys = dense.(ys); dzs = dense.(zs)
    push!(rs, run_per_op(name, M, () -> map(op, xs, ys, zs), () -> map(op, dxs, dys, dzs),
        (c, r) -> densematch(c, r; rtol = RT)))
end

#---KVector workloads: geometric primitives & incidence--------------------------------------------#

function kvector_benchmarks()
    rng = MersenneTwister(0x6EC7_0123)
    M = 1000
    rs = WResult[]

    # PGA(3): point = grade 3, line = grade 2, plane = grade 1.
    planes = [randk(rng, 1, PGA(3)) for _ in 1:M]
    planes2 = [randk(rng, 1, PGA(3)) for _ in 1:M]
    planes3 = [randk(rng, 1, PGA(3)) for _ in 1:M]
    lines = [randk(rng, 2, PGA(3)) for _ in 1:M]
    points = [randk(rng, 3, PGA(3)) for _ in 1:M]
    points2 = [randk(rng, 3, PGA(3)) for _ in 1:M]

    cvd2!(rs, "meet: plane ∧ plane → line (PGA3)", M, wedge, planes, planes2)
    cvd2!(rs, "meet: line ∧ plane → point (PGA3)", M, wedge, lines, planes)
    cvd3!(rs, "meet: 3 planes → point (PGA3)", M, (a, b, c) -> a ∧ b ∧ c, planes, planes2, planes3)
    cvd2!(rs, "join: point ∨ point → line (PGA3)", M, regressive, points, points2)
    cvd2!(rs, "join: line ∨ point → plane (PGA3)", M, regressive, lines, points)
    cvd2!(rs, "project: (π⌋P)∧π → point (PGA3)", M,
        (pl, p) -> wedge(left_contraction(pl, p), pl), planes, points)

    # PGA(2): point = grade 2, line = grade 1.
    plines = [randk(rng, 1, PGA(2)) for _ in 1:M]
    plines2 = [randk(rng, 1, PGA(2)) for _ in 1:M]
    ppoints = [randk(rng, 2, PGA(2)) for _ in 1:M]
    ppoints2 = [randk(rng, 2, PGA(2)) for _ in 1:M]
    cvd2!(rs, "meet: line ∧ line → point (PGA2)", M, wedge, plines, plines2)
    cvd2!(rs, "join: point ∨ point → line (PGA2)", M, regressive, ppoints, ppoints2)

    # VGA(3): direction = grade 1, plane element / rotation generator = grade 2.
    vs = [randk(rng, 1, VGA(3)) for _ in 1:M]
    ws = [randk(rng, 1, VGA(3)) for _ in 1:M]
    bvs = [randk(rng, 2, VGA(3)) for _ in 1:M]
    bvs2 = [randk(rng, 2, VGA(3)) for _ in 1:M]
    cvd2!(rs, "wedge: v ∧ w → bivector (VGA3)", M, wedge, vs, ws)
    cvd2!(rs, "scalar_product ⟨v,w⟩ (VGA3)", M, scalar_product, vs, ws)
    cvd2!(rs, "left_contraction v⌋B (VGA3)", M, left_contraction, vs, bvs)
    cvd1!(rs, "right_complement bivector (VGA3)", M, right_complement, bvs)
    cvd2!(rs, "commutator B₁×B₂ (VGA3)", M, commutator, bvs, bvs2)
    cvd1!(rs, "reverse line bivector (PGA3)", M, reverse, lines)

    return rs
end

#---EvenCliffordNumber workloads: rotors & motors--------------------------------------------------#

function even_benchmarks()
    rng = MersenneTwister(0x4EE7_0456)
    M = 1000
    rs = WResult[]

    # VGA(3) rotors acting on grade-1 vectors.
    r3 = [randversor(rng, VGA(3)) for _ in 1:M]
    r3b = [randversor(rng, VGA(3)) for _ in 1:M]
    bvs = [randk(rng, 2, VGA(3)) for _ in 1:M]
    vs = [randk(rng, 1, VGA(3)) for _ in 1:M]
    cvd2!(rs, "rotor compose r₁r₂ (VGA3)", M, *, r3, r3b)
    cvd2!(rs, "rotor sandwich r·v·r̃ (VGA3)", M, sandwich, r3, vs)
    cvd1!(rs, "rotor from bivector exp(B) (VGA3)", M, exp, bvs)
    cvd2!(rs, "relative rotor r₁/r₂ (VGA3)", M, /, r3, r3b)

    # PGA(3) motors (rigid motions) acting on points / lines / planes.
    motors = [randversor(rng, PGA(3)) for _ in 1:M]
    motors2 = [randversor(rng, PGA(3)) for _ in 1:M]
    pbvs = [randk(rng, 2, PGA(3)) for _ in 1:M]
    pts = [randk(rng, 3, PGA(3)) for _ in 1:M]
    lns = [randk(rng, 2, PGA(3)) for _ in 1:M]
    pls = [randk(rng, 1, PGA(3)) for _ in 1:M]
    cvd2!(rs, "motor compose m₁m₂ (PGA3)", M, *, motors, motors2)
    cvd2!(rs, "move point M·P·M̃ (PGA3)", M, sandwich, motors, pts)
    cvd2!(rs, "move line M·L·M̃ (PGA3)", M, sandwich, motors, lns)
    cvd2!(rs, "move plane M·π·M̃ (PGA3)", M, sandwich, motors, pls)
    cvd1!(rs, "motor from bivector exp(B) (PGA3)", M, exp, pbvs)

    # CGA(3) motor acting on a grade-1 conformal point.
    cmotors = [randeven2(rng, CGA(3)) for _ in 1:M]
    cpts = [randk(rng, 1, CGA(3)) for _ in 1:M]
    cvd2!(rs, "CGA move point M·P·M̃ (CGA3)", M, sandwich, cmotors, cpts)

    # VGA(4) rotor (so(4) double rotation) composition.
    r4 = [randversor(rng, VGA(4)) for _ in 1:M]
    r4b = [randversor(rng, VGA(4)) for _ in 1:M]
    cvd2!(rs, "rotor compose (VGA4)", M, *, r4, r4b)

    return rs
end

#---OddCliffordNumber workloads: reflections------------------------------------------------------#

function odd_benchmarks()
    rng = MersenneTwister(0x0DD0_0789)
    M = 1000
    rs = WResult[]

    # Grade-1 versor (plane / vector) reflections: the odd-grade sandwich.
    planes = [randk(rng, 1, PGA(3)) for _ in 1:M]
    planes2 = [randk(rng, 1, PGA(3)) for _ in 1:M]
    points = [randk(rng, 3, PGA(3)) for _ in 1:M]
    cvd2!(rs, "reflect point π·P·π (PGA3)", M, (pl, p) -> pl * p * pl, planes, points)
    cvd2!(rs, "reflect plane π·π′·π (PGA3)", M, (pl, q) -> pl * q * pl, planes, planes2)

    vs = [randk(rng, 1, VGA(3)) for _ in 1:M]
    ns = [randk(rng, 1, VGA(3)) for _ in 1:M]
    cvd2!(rs, "reflect vector n·v·n (VGA3)", M, (n, v) -> n * v * n, ns, vs)

    # Genuine OddCliffordNumber versors (product of three 1-vectors).
    odds = [randodd(rng, VGA(3)) for _ in 1:M]
    odds2 = [randodd(rng, VGA(3)) for _ in 1:M]
    rotors = [randversor(rng, VGA(3)) for _ in 1:M]
    vs2 = [randk(rng, 1, VGA(3)) for _ in 1:M]
    cvd2!(rs, "odd versor sandwich u·v·u⁻¹ (VGA3)", M,
        (u, v) -> u * v * versor_inverse(u), odds, vs2)
    cvd2!(rs, "odd∘odd → rotor u₂u₁ (VGA3)", M, *, odds, odds2)
    cvd2!(rs, "rotor∘reflection r·u → odd (VGA3)", M, *, rotors, odds)
    cvd1!(rs, "grade_involution odd versor (VGA3)", M, grade_involution, odds)

    return rs
end

#---General CliffordNumber workloads: full multivectors & inverses---------------------------------#

function general_benchmarks()
    rng = MersenneTwister(0x9E11_0ABC)
    M = 1000
    rs = WResult[]

    # Full geometric product in VGA(4) timed against an even/odd block decomposition
    # of the same product: an algebraic identity (x = x₊ + x₋), comparing the full
    # product against the summed graded blocks.
    xs = [CliffordNumber{VGA(4),Float64}(ntuple(_ -> randn(rng), 16)) for _ in 1:M]
    ys = [CliffordNumber{VGA(4),Float64}(ntuple(_ -> randn(rng), 16)) for _ in 1:M]
    xes = EvenCliffordNumber{VGA(4)}.(xs); xos = OddCliffordNumber{VGA(4)}.(xs)
    yes = EvenCliffordNumber{VGA(4)}.(ys); yos = OddCliffordNumber{VGA(4)}.(ys)
    blockmul(xe, xo, ye, yo) = dense(xe * ye) + dense(xe * yo) + dense(xo * ye) + dense(xo * yo)
    push!(rs, run_per_op("geometric product full vs block (VGA4)", M,
        () -> map(*, xs, ys),
        () -> map(blockmul, xes, xos, yes, yos),
        (c, r) -> approxall(c, r; rtol = RT)))

    # `inv` validates `x·inv ≈ 1`; `versor_inverse` skips the check. For a true versor
    # the results are equal, so the ratio is the cost of that validation, here for
    # rotors and motors.
    r3 = [randversor(rng, VGA(3)) for _ in 1:M]
    motors = [randversor(rng, PGA(3)) for _ in 1:M]
    push!(rs, run_per_op("inv vs versor_inverse: rotor (VGA3)", M,
        () -> map(inv, r3), () -> map(versor_inverse, r3), (c, r) -> approxall(c, r; rtol = RT)))
    push!(rs, run_per_op("inv vs versor_inverse: motor (PGA3)", M,
        () -> map(inv, motors), () -> map(versor_inverse, motors), (c, r) -> approxall(c, r; rtol = RT)))

    return rs
end

#---Batched-broadcast SIMD workloads--------------------------------------------------------------#
#
# Apply one fixed transform / operation across an array of primitives at increasing N
# (the vertex-buffer / point-cloud workload). Compact vs dense, with the drift check
# sampled over the leading slice so the 1M-element rows don't allocate a second dense
# copy of the whole array just to verify. A flat compact/dense ratio across N means
# the kernel scales; a ratio that grows with N flags a per-element path that fails to
# vectorize or turns memory-bound.

const SIMD_SIZES = (1_024, 65_536, 1_048_576)

# Sampled compact-vs-dense check: all elements run the same kernel, so a leading
# slice is a sufficient drift oracle without densifying a million-element array.
function sampled_match(cs, rs; n = 256)
    k = min(n, length(cs), length(rs))
    return densematch(view(cs, 1:k), view(rs, 1:k); rtol = RT)
end

# A fixed transform `R` broadcast over an array of primitives via the sandwich `R X R̃`.
function simd_sandwich!(rs, label, R, xs)
    dR = dense(R); dxs = dense.(xs)
    for N in SIMD_SIZES
        v = @view xs[1:N]; dv = @view dxs[1:N]
        push!(rs, run_workload("$label (N=$N)",
            () -> sandwich.(Ref(R), v),
            () -> sandwich.(Ref(dR), dv),
            (c, r) -> sampled_match(c, r)))
    end
end

# A binary op broadcast element-wise over two arrays.
function simd_binary!(rs, label, op, xs, ys)
    dxs = dense.(xs); dys = dense.(ys)
    for N in SIMD_SIZES
        vx = @view xs[1:N]; vy = @view ys[1:N]
        dvx = @view dxs[1:N]; dvy = @view dys[1:N]
        push!(rs, run_workload("$label (N=$N)",
            () -> op.(vx, vy),
            () -> op.(dvx, dvy),
            (c, r) -> sampled_match(c, r)))
    end
end

function simd_benchmarks()
    rng = MersenneTwister(0x5114_0DEF)
    Nmax = maximum(SIMD_SIZES)
    rs = WResult[]

    # Rigid-motion of a point cloud / vertex buffer (PGA3).
    M = randversor(rng, PGA(3))
    pts = [randk(rng, 3, PGA(3)) for _ in 1:Nmax]
    simd_sandwich!(rs, "move point M·P·M̃ (PGA3)", M, pts)

    # Reflection of a point cloud across a fixed plane (PGA3): odd-grade sandwich.
    plane = randk(rng, 1, PGA(3))
    pts2 = [randk(rng, 3, PGA(3)) for _ in 1:Nmax]
    simd_sandwich!(rs, "reflect point π·P·π (PGA3)", plane, pts2)

    # Rotation of a direction array (VGA3).
    r = randversor(rng, VGA(3))
    vs = [randk(rng, 1, VGA(3)) for _ in 1:Nmax]
    simd_sandwich!(rs, "rotor sandwich r·v·r̃ (VGA3)", r, vs)

    # Incidence over arrays: meet of plane pairs and join of point pairs (PGA3).
    pa = [randk(rng, 1, PGA(3)) for _ in 1:Nmax]
    pb = [randk(rng, 1, PGA(3)) for _ in 1:Nmax]
    simd_binary!(rs, "meet plane ∧ plane (PGA3)", wedge, pa, pb)
    P1 = [randk(rng, 3, PGA(3)) for _ in 1:Nmax]
    P2 = [randk(rng, 3, PGA(3)) for _ in 1:Nmax]
    simd_binary!(rs, "join point ∨ point (PGA3)", regressive, P1, P2)

    # Rotor composition over arrays (VGA3): pure even×even geometric product.
    ra = [randversor(rng, VGA(3)) for _ in 1:Nmax]
    rb = [randversor(rng, VGA(3)) for _ in 1:Nmax]
    simd_binary!(rs, "rotor compose r₁r₂ (VGA3)", *, ra, rb)

    return rs
end
