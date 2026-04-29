# Structured Poisson: LDU vs DIA vs RB-DIA Gauss-Seidel Benchmark

A standalone C/C++ benchmark comparing three storage/ordering schemes for the
Gauss-Seidel smoother on a 3D Poisson problem discretised on a structured
hexahedral mesh:

- **LDU** — indirect addressing via (lower, diag, upper) face lists, mirroring
  OpenFOAM's native `lduMatrix` layout
- **DIA** — direct diagonal storage, using the fixed offset stencil {1, N, N²}
  of a structured hex mesh to eliminate indirect addressing
- **RB-DIA** — Red-Black ordered DIA variant with OpenMP parallelism across
  cells of the same colour, enabling shared-memory multi-core execution

The problem is a 3D manufactured-solution Poisson equation
`-∇²u = 3π² sin(πx) sin(πy) sin(πz)` with homogeneous Dirichlet BCs, discretised
on an N×N×N unit cube. All three solvers produce bit-for-bit identical
converged solutions, verified by matching L2 error against the manufactured
exact solution with second-order mesh convergence.

## Why this benchmark exists

This repo is the standalone testbed backing a larger OpenFOAM plugin project
([DIAGaussSeidel-Smoother-OpenFOAM](https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM)),
which replaces OpenFOAM's stock Gauss-Seidel smoother with a DIA-based variant
for structured hex meshes. Inside OpenFOAM the plugin yields ~50% reduction in
cache misses at all levels and a ~28% speedup on the GAMG smoother kernel.

The benchmark here serves two purposes:

1. **Clean, reproducible measurement of the core kernel** outside the OpenFOAM
   framework, useful for profiling under LIKWID and for studying algorithmic
   properties without the surrounding simulation machinery.

2. **A portable platform** for running the same kernel on non-x86 hardware —
   specifically Calligo Technologies' **Uttunga** RISC-V processor with native
   posit arithmetic. On Uttunga, the goal shifts from cache efficiency (which
   motivated the OpenFOAM plugin) to studying whether the higher precision of
   posit arithmetic reduces iteration counts to convergence — analogous to
   results previously observed with HYPRE BoomerAMG on multimaterial thermal
   analysis (aluminium + CFRP).

## Modes

The benchmark uses an interactive menu. Run the binary and select a mode:

| Mode | Menu option | Purpose |
|---|---|---|
| Kernel profiling | 1 | Times a fixed number of sweeps for LDU, DIA, and RBDIA with optional LIKWID markers. Default: N=50, 1000 sweeps |
| Combined convergence | 2 | Runs DIA then RBDIA convergence tests back-to-back. Default: N=50, tol=1e-16, max_sweeps=50000 |
| DIA convergence | 3 | Iterates DIA until residual drops below tolerance, reports sweep count and L2 error |
| RBDIA convergence | 4 | Iterates RBDIA until residual drops below tolerance, reports sweep count and L2 error. Set `OMP_NUM_THREADS` for parallelism |
| Math sanity check | 5 | Validates math library functions, grid coordinate computation, and source assembly across platforms |
| Help | 6 | Prints usage information and defaults |

All prompts accept Enter to use the default value (N=50, tol=1e-16, max_sweeps=50000).

For unattended runs (e.g. over SSH to Uttunga), pipe inputs:
```
echo -e "3\n50\n1e-16\n50000" | ./poisson_benchmark 2>&1 | tee results.log
```

## Build

Requires a C++14 compiler and OpenMP.

```
make
```

For LIKWID profiling, compile with:
```
make profile
```
This enables `-DLIKWID_PERFMON` and links against `-llikwid`.

## Results

### x86 kernel performance (AMD Zen4, 8 cores)

| Solver | 1000-sweep time at N=100 (s) | OMP scaling (8 threads) |
|---|---|---|
| DIA (serial) | ~7.5 | n/a |
| RB-DIA (1 thread) | ~8.2 | baseline |
| RB-DIA (2 threads) | ~1.7 | 1.82× (91% efficiency) |
| RB-DIA (4 threads) | ~1.0 | 2.95× (74% efficiency) |
| RB-DIA (8 threads) | ~2.1 | 3.89× (49% efficiency) |

DIA is inherently serial due to the forward-push data dependency in
natural-order Gauss-Seidel. RB-DIA breaks that dependency via Red-Black
ordering and scales well under OpenMP, with near-linear speedup at low
thread counts.

### Posit vs IEEE 754: convergence comparison (N=50, tol=1e-16)

| Solver | Platform | Sweeps to converge | Final residual | L2 error |
|---|---|---|---|---|
| DIA | Uttunga (posit) | **7,950** | 9.67e-17 | 1.163374e-04 |
| DIA | x86 (IEEE 754) | **50,000+ (failed)** | 1.95e-16 (stalled) | 1.163374e-04 |
| RB-DIA | Uttunga (posit) | **8,050** | 9.18e-17 | 1.163374e-04 |
| RB-DIA | x86 (IEEE 754) | **50,000+ (failed)** | 1.97e-16 (stalled) | 1.163374e-04 |

On x86, IEEE 754 double precision limits the achievable residual to ~2e-16
(the machine epsilon floor). The solver stalls at this level regardless of
additional iterations. On Uttunga with posit arithmetic, the solver drives
the residual below this floor and converges in ~8,000 sweeps.

Both platforms produce identical L2 errors (1.163374e-04), confirming they
converge to the same discrete solution. The difference is purely in how
precisely each number system can solve the algebraic system Ax = b.

This result is consistent with prior observations on HYPRE BoomerAMG, where
posit arithmetic on Uttunga reduced iteration counts on ill-conditioned
multimaterial thermal problems (aluminium + CFRP). The present benchmark
reproduces the same phenomenon on a simpler, well-conditioned Poisson
problem using a different solver class (stationary Gauss-Seidel vs
algebraic multigrid).

### Hardware

- **x86:** AMD Ryzen 7 250 (Zen4), 8 cores, 32 MB L3, Ubuntu 24.04, GCC 11, `-O3`
- **Uttunga:** Calligo Technologies Uttunga v1.0, RISC-V with posit arithmetic
  co-processor, 4 cores at 400 MHz, 2 GB RAM

## Related work

- [DIAGaussSeidel-Smoother-OpenFOAM](https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM) —
  the OpenFOAM Foundation v13 plugin that deploys the DIA storage scheme as a
  drop-in GAMG smoother, validated on turbulent cavity flows with ~50% cache
  miss reduction and ~28% kernel speedup
- [RBDIAGaussSeidel-Smoother-OpenFOAM](https://github.com/amartyadav/RBDIAGaussSeidel-Smoother-OpenFOAM) —
  the Red-Black DIA variant with OpenMP, designed for offload to Uttunga's
  multi-core architecture via Calligo's PETSc-OpenFOAM integration

## Author

Amartya Yadav <br/>
HPC Software Engineer,<br/> Calligo Technologies

## License

MIT
