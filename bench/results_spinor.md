# PR-SpinorFastPaths — before/after

Effect of the Group B fast paths on the [`bench/`](README.md) suite:

1. **fast `inv` / `versor_inverse`** for the ℝ/ℂ/ℍ even subalgebras — a single
   fused closed-form reciprocal (skip the validating product, fold in the reverse,
   reciprocate the modulus once instead of dividing each component),
2. **closed-form `*`/`²`** for the 2-component even spinor, and
3. **closed-form `abs2`** for the positive-definite even subalgebras.

**Before** is `main` (the recorded run in [`results.md`](results.md)); **after**
is branch `pr-spinorfastpaths` with all three. Both runs used `BENCH_SECONDS=5`
on the same machine, and every workload still passes its algebra-drift check —
the fast paths are pure speedups, not behavior changes.

`CN/ref` is the CliffordNumbers time over the matching Julia primitive
(`Complex{T}` / `Quaternion{T}`); lower is better, `<1` means CliffordNumbers is
faster. The bench suite itself is the soft dependency PR-BenchSuite (not part of
this branch); these numbers were recorded against a local copy.

## Environment

| | |
|---|---|
| Branch | `pr-spinorfastpaths` (baseline `main` @ `c5ca350`) |
| Date | 2026-06-04 |
| Julia | 1.12.5 |
| CPU | x86_64 skylake, 16 threads |
| OS | Windows 11 (`x86_64-w64-mingw32`) |
| `BENCH_SECONDS` | 5.0 |

## bench_complex — even VGA(2) ≅ ℂ

| workload | before CN | before CN/ref | after CN | after CN/ref | speedup |
|---|---:|---:|---:|---:|---:|
| **microbench: `inv` (per op)** | 56.2 ns | 20.81 | **1.1 ns** | **0.41** | **51×** |
| **contour integral (N=4096)** | 675.3 µs | 17.63 | **18.3 µs** | **0.48** | **37×** |
| **Mandelbrot 512² (1 thread)** | 18.359 ms | 2.03 | **8.603 ms** | **0.97** | **2.1×** |
| **Mandelbrot 512² (threaded)** | 2.717 ms | 2.14 | **1.214 ms** | **0.97** | **2.2×** |
| **Julia set 512²** | 13.899 ms | 2.10 | **6.222 ms** | **0.95** | **2.2×** |
| **microbench: `versor_inverse` (per op)** | 9.8 ns | 3.63 | **1.1 ns** | **0.42** | **8.9×** |
| microbench: `*` (per op) | 0.9 ns | 1.00 | 0.9 ns | 1.00 | 1× |
| pointwise `*` (N=65536) | 160.5 µs | 1.00 | 162.9 µs | 1.03 | ≈1× |
| `z² + c` (N=65536) | 161.5 µs | 1.00 | 158.4 µs | 0.97 | ≈1× |
| Horner deg-8 (N=65536) | 590.4 µs | 0.90 | 594.2 µs | 0.95 | ≈1× |

`inv` and `versor_inverse` share the closed form (`inv` delegates to it), so both
now beat the primitive at ≈0.4×.

### Closed-form reciprocal, focused (suite has no ℍ-`inv` row)

Per-op microbench (`map(inv, xs)` over 1000 elements), comparing the closed-form
`inv` / `versor_inverse` to the previous path and the primitive:

| `inv` per op | previous (generic versor) | closed form | primitive ref |
|---|---:|---:|---:|
| spinor (VGA(2) ≅ ℂ) | 4.4 ns | **1.1 ns** | `Complex.inv` 2.7 ns |
| rotor (VGA(3) ≅ ℍ) | 11.1 ns | **2.1 ns** | `Quaternion.inv` 29.2 ns |

The ℍ case is where collapsing **4 component divides → 1 reciprocal + 4 muls**
pays off most: ~5× over the old path, ~14× over `Quaternion.inv`. Both Base and
`Quaternions.jl` `inv` do overflow-scaling that the versor formula skips, which
is why the closed form lands below the primitive for normal-range inputs.

## bench_quaternion — even VGA(3) ≅ ℍ

The ℍ workloads here exercise the geometric product (`compose`, `sandwich`,
`chain`), not `inv`/`abs2`, and the closed-form `*` covers only the 2-component
spinor — so they are unchanged. Shown for completeness; differences are noise.

| workload | before CN/ref | after CN/ref |
|---|---:|---:|
| compose (per op) | 2.89 | 2.89 |
| sandwich (per op) | 0.74 | 0.73 |
| batched sandwich (N=100) | 0.74 | 0.74 |
| batched sandwich (N=100k) | 0.75 | 0.76 |
| batched sandwich (N=1M) | 0.79 | 0.79 |
| slerp (per op) | 1.04 | 1.01 |
| chain-compose (N=10k) | 0.75 | 0.77 |

## What moved, and why

The three changes hit three different bottlenecks; together the spinor is now
**at parity-or-better than `Complex{T}` on every ℂ workload**, and `inv` /
contour are well below it.

**Closed-form reciprocal.** Beyond dropping the validating `x·inv ≈ 1` check,
emitting the inverse as one fused expression — reverse folded into the signs (no
intermediate `x'`), `abs2` inlined, modulus reciprocated once and multiplied
through (1 divide vs 2 / 4) — takes `inv` from **20.81× to 0.41×** and the contour
integral (its integrand is `1/(z − pₖ)`) from **17.63× to 0.48×**. The closed form
lives in `versor_inverse` (so the public versor inverse is fast too,
**3.63× → 0.42×**), and `inv` delegates to it.

**Closed-form `abs2`.** `abs2(x) = scalar_product(x, x')` built the reverse and
ran a three-tuple `mapreduce`; for the positive-definite even subalgebras every
blade has `sign_of_square·reverse_sign = +1`, so it is a plain sum of squares.
This was the dominant win for the escape loops — the `abs2(z) > 4` test runs
every iteration, so Mandelbrot/Julia went **~2.1× → ~0.97×** (now faster than
`Complex{T}`).

**Closed-form `*`/`²`.** Removes the un-foldable `muladd(0,0,x)` the generic
kernel seeds with. It matters only where the multiply can't vectorize: the scalar
escape loops improve, while the already-vectorized array rows (`pointwise *`,
`z² + c`) sit at parity throughout.

## Remaining hot spots

- **ℍ compose stays ~2.8×, and a closed-form Hamilton product does not help.** It
  was tried (both a scalar `muladd` chain and a SIMD-friendly broadcast form) and
  dropped: unlike the 2-component spinor, the generic `@generated` kernel already
  SIMD-vectorizes the 4-component product over 4-tuples, so removing its
  `zero_tuple` seed gave nothing, and the closed forms measured *slower* in the
  compute-bound reduce case (`chain-compose` ≈ 5.3 ns/mul generic vs ≈ 7.4–7.9
  ns/mul closed). The ~2.8× gap to `Quaternion` is not the geometric product
  itself — `Quaternions.jl` `*` is a hand-tuned 4-term form, while the GA value
  passes through promotion/`muladd` machinery that is already well vectorized.

## Test coverage

The specialized methods (`inv`, `*`, `abs2` on the even subalgebras) are covered
by the `@testset "Even subalgebra fast paths"` in `test/operations.jl`: exact
integer cases pin the product formula across all three σ branches (ℂ / split /
dual), and randomized cases cross-check the closed forms against the generic
kernel, `versor_inverse`, `scalar_product`, and `Complex`.

## Reproducing

```sh
# against this branch, with the PR-BenchSuite harness in bench/
BENCH_SECONDS=5 julia -t auto bench/runbenchmarks.jl
```
