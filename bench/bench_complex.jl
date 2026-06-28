# Spinor (≅ ℂ) workloads — the acceptance baseline for PR-SpinorFastPaths.
#
# `Spinor{T} = EvenCliffordNumber{VGA(2),T,2}` is isomorphic to `Complex{T}`
# coordinate-for-coordinate, so every workload compares the spinor result to the
# `Complex{T}` result with `≈` (drift check) and times one against the other.

# Generic kernels — these compile identically for `Spinor` and `Complex`.

function horner(x, c::NTuple{N}) where N
    acc = c[N] * one(x)
    @inbounds for k in (N - 1):-1:1
        acc = acc * x + c[k]
    end
    return acc
end

escape_count(z0, c, maxiter) = begin
    z = z0
    @inbounds for i in 1:maxiter
        z = z^2 + c
        abs2(z) > 4 && return i
    end
    return maxiter
end

function escape_grid(point, xs, ys, z0_of, c_of, maxiter; threaded::Bool=false)
    n = length(xs)
    m = length(ys)
    out = Matrix{Int}(undef, n, m)
    if threaded
        Threads.@threads for j in 1:m
            @inbounds for i in 1:n
                p = point(xs[i], ys[j])
                out[i, j] = escape_count(z0_of(p), c_of(p), maxiter)
            end
        end
    else
        @inbounds for j in 1:m, i in 1:n
            p = point(xs[i], ys[j])
            out[i, j] = escape_count(z0_of(p), c_of(p), maxiter)
        end
    end
    return out
end

# Trapezoidal contour integral of Σ 1/(z - pₖ) around a closed loop of samples.
# The `dz` factor makes the result ≈ 2πi·(enclosed poles) — a stable nonzero
# reference, unlike the bare residue sum, which cancels to ~0 on a symmetric loop.
function contour_integral(pts, poles)
    n = length(pts)
    s = zero(pts[1])
    @inbounds for i in 1:n
        z = pts[i]
        dz = pts[i == n ? 1 : i + 1] - z
        f = zero(pts[1])
        for p in poles
            f += inv(z - p)
        end
        s += f * dz
    end
    return s
end

#---Workload assembly------------------------------------------------------------------------------#

function complex_benchmarks()
    rng = MersenneTwister(0xC1FF)
    T = Float64
    N = 65_536
    G = 512          # Mandelbrot / Julia grid edge
    MAXITER = 32

    results = WResult[]
    approxall(a, b) = all(((x, y),) -> isapprox(x, y; rtol = sqrt(eps(T))), zip(a, b))

    # --- pointwise *, z^2 + c, Horner over N points ---
    zs = [complex(2randn(rng) - 1, 2randn(rng) - 1) for _ in 1:N]
    ws = [complex(2randn(rng) - 1, 2randn(rng) - 1) for _ in 1:N]
    cs = [complex(0.3randn(rng), 0.3randn(rng)) for _ in 1:N]
    szs = tospinor.(zs); sws = tospinor.(ws); scs = tospinor.(cs)
    coeffs = ntuple(_ -> randn(rng), 9)   # degree-8 polynomial

    push!(results, run_workload("pointwise *  (N=65536)",
        () -> szs .* sws,
        () -> zs .* ws,
        (c, r) -> approxall(ascomplex.(c), r)))

    push!(results, run_workload("z^2 + c  (N=65536)",
        () -> szs .^ 2 .+ scs,
        () -> zs .^ 2 .+ cs,
        (c, r) -> approxall(ascomplex.(c), r)))

    push!(results, run_workload("Horner deg-8  (N=65536)",
        () -> horner.(szs, Ref(coeffs)),
        () -> horner.(zs, Ref(coeffs)),
        (c, r) -> approxall(ascomplex.(c), r)))

    # --- Mandelbrot 512^2, single-threaded and threaded ---
    mxs = range(-2.0, 0.5; length = G)
    mys = range(-1.25, 1.25; length = G)
    spt(re, im) = Spinor{T}((re, im))
    cpt(re, im) = complex(re, im)
    push!(results, run_workload("Mandelbrot 512^2 (1 thread)",
        () -> escape_grid(spt, mxs, mys, zero, identity, MAXITER),
        () -> escape_grid(cpt, mxs, mys, zero, identity, MAXITER),
        (c, r) -> c == r))

    push!(results, run_workload("Mandelbrot 512^2 (threaded)",
        () -> escape_grid(spt, mxs, mys, zero, identity, MAXITER; threaded = true),
        () -> escape_grid(cpt, mxs, mys, zero, identity, MAXITER; threaded = true),
        (c, r) -> c == r))

    # --- Julia set 512^2 (fixed c, varying seed point) ---
    jxs = range(-1.5, 1.5; length = G)
    jys = range(-1.5, 1.5; length = G)
    sjc = Spinor{T}((-0.8, 0.156))
    cjc = complex(-0.8, 0.156)
    push!(results, run_workload("Julia set 512^2",
        () -> escape_grid(spt, jxs, jys, identity, _ -> sjc, MAXITER),
        () -> escape_grid(cpt, jxs, jys, identity, _ -> cjc, MAXITER),
        (c, r) -> c == r))

    # --- contour integral with poles, N sample points on a circle ---
    NC = 4096
    θs = range(0, 2π; length = NC + 1)[1:NC]
    ring = [complex(3cos(θ), 3sin(θ)) for θ in θs]
    poles = (complex(0.5, 0.5), complex(-1.0, 0.3), complex(0.2, -1.4))
    sring = tospinor.(ring); spoles = tospinor.(poles)
    push!(results, run_workload("contour integral (N=4096)",
        () -> contour_integral(sring, spoles),
        () -> contour_integral(ring, poles),
        (c, r) -> isapprox(ascomplex(c), r; rtol = sqrt(eps(T)))))

    # --- per-op microbenches: the per-call acceptance baselines ---
    # Batched over a small array so the op can't be constant-folded; reported
    # times are per element (see `run_per_op`).
    M = 1000
    mz = [complex(2randn(rng) - 1, 2randn(rng) - 1) for _ in 1:M]
    ms = tospinor.(mz)
    sq(x) = x * x
    push!(results, run_per_op("microbench: * (per op)", M,
        () -> map(sq, ms), () -> map(sq, mz), (c, r) -> approxall(ascomplex.(c), r)))
    push!(results, run_per_op("microbench: inv (per op)", M,
        () -> map(inv, ms), () -> map(inv, mz), (c, r) -> approxall(ascomplex.(c), r)))
    push!(results, run_per_op("microbench: versor_inverse (per op)", M,
        () -> map(versor_inverse, ms), () -> map(inv, mz), (c, r) -> approxall(ascomplex.(c), r)))

    return results
end
