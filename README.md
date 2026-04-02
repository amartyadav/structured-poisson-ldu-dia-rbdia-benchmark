# Structured Gauss-Seidel Benchmark: LDU vs DIA Addressing

A standalone C benchmark comparing OpenFOAM's indirect-addressed (LDU) Gauss-Seidel smoother against a diagonal-format (DIA) implementation with compile-time-known offsets, targeting structured hexahedral meshes. This is **Phase 1** of a larger project to optimise OpenFOAM's linear solver stack for structured meshes.

## Background and Motivation

OpenFOAM's `lduMatrix` format stores sparse matrices using face-based indirect addressing: two index arrays (`upperAddr`, `lowerAddr`) map each internal face to the pair of cells it connects. The off-diagonal coefficients (`upper[]`, `lower[]`) are accessed via these indirection arrays during every Gauss-Seidel sweep, SpMV (`Amul`), and preconditioner application.

This design is fully general — it handles any mesh topology. But on structured hex meshes (the workhorse of industrial CFD), every cell has the same neighbour offsets. A 50×50×50 cavity mesh has exactly three unique offsets: `{±1, ±50, ±2500}`. The indirection arrays are redundant: they encode a pattern that could be expressed as compile-time constants.

The cost of this redundancy is severe. Phase 0 profiling (LIKWID, AMD Ryzen 7 250/Zen4) on a 50³ cavity case showed:

- **0.31% of peak FLOP/s** (1.93 GFLOP/s out of ~627 GFLOP/s theoretical) — cores are starved for data
- **7–8% L1D miss rate** on the hot threads running the smoother
- **38–50% L2 miss rate** — half of everything that misses L1 also misses L2
- **46.5% L3 miss rate** — nearly half of L3 lookups go all the way to DRAM
- Estimated **~9 GB/s DRAM bandwidth** out of ~50–70 GB/s peak (13–18% utilisation)

The root cause is the indirect `psi[uAddr[f]]` access pattern in the GS sweep inner loop. The hardware prefetcher cannot predict the next address because it depends on a value loaded from the index array. Every face access is effectively a pointer chase.

The DIA (Diagonal) format eliminates this entirely. Instead of `upper[]` indexed by faces with indirect `uAddr[]` lookups, the matrix is stored as contiguous coefficient arrays aligned to each diagonal offset. The GS sweep becomes `psi[i±1]`, `psi[i±N]`, `psi[i±N²]` — fully strided, prefetchable, and vectorisable.

This benchmark isolates the addressing scheme as the single variable. Both implementations solve the identical problem, use the identical algorithm (forward Gauss-Seidel), and produce the identical answer. The only difference is how off-diagonal coefficients and neighbour values are looked up.

## The Larger Project

This repo is one phase of a three-phase project at Calligo Technologies:

| Phase | Description | Status |
|-------|-------------|--------|
| **Phase 0** | LIKWID profiling baseline on OpenFOAM 50³ cavity; lduMatrix layout study; diagonal offset verification; AMD IBS analysis of `renumberMesh` effects | Complete |
| **Phase 1** | Standalone benchmark (this repo) — head-to-head LDU vs DIA Gauss-Seidel on a 2D Poisson problem with manufactured solution | In progress |
| **Phase 2** | OpenFOAM drop-in smoother plugin — `structuredGaussSeidelSmoother` registered via runtime selection, with structured mesh auto-detection and LDU→DIA conversion at setup time | Planned |

The kernel developed here is designed to be portable for a future port to the Uttunga RISC-V posit co-processor board.

### Key Phase 0 Findings

Phase 0 established two critical facts beyond the profiling numbers above:

1. **Fixed diagonal offsets are real.** Instrumentation of `GaussSeidelSmoother.C` confirmed that for the 50³ structured hex cavity, the only neighbour offsets produced by `lduAddressing` are `{1, 50, 2500}` — no exceptions. The DIA reformulation is exact, not an approximation.

2. **Cuthill-McKee reordering helps SpMV but hurts the smoother.** AMD IBS profiling showed that `renumberMesh` (which applies Cuthill-McKee bandwidth reduction) improved `Foam::Amul` performance but degraded `Foam::smooth`. The reason: CMK reordering makes `psi[uPtr[facei]]` indices less predictable relative to `celli`, increasing cache miss scatter in the GS sweep even as it improves sequential face traversal for SpMV. The DIA format sidesteps this entirely — offsets are fixed constants regardless of cell numbering.

## Problem Formulation

The benchmark solves the 2D Poisson equation on the unit square with homogeneous Dirichlet boundary conditions:

```
−(∂²u/∂x² + ∂²u/∂y²) = f(x,y)     on [0,1] × [0,1]
                    u = 0            on all boundaries
```

The source term is chosen to produce a known exact solution (method of manufactured solutions):

```
f(x,y) = 2π² · sin(πx) · sin(πy)
u_exact(x,y) = sin(πx) · sin(πy)
```

Discretisation uses second-order central differences on a uniform N×N cell-centred grid, producing the standard 5-point stencil. Cell `(i,j)` maps to flat index `i*N + j` (row-major). The resulting linear system `Ax = b` has `N²` unknowns with at most 4 off-diagonal entries per row (fewer for boundary cells).

## Implementation Details

### LDU Solver (`ldu_solver.c`)

Mirrors OpenFOAM's `GaussSeidelSmoother::smooth()` exactly:

- **Storage:** Five flat arrays — `diag[nCells]`, `upper[nFaces]`, `lower[nFaces]`, `upperAddr[nFaces]`, `lowerAddr[nFaces]` — plus an `ownerStart[nCells+1]` array (CSR-style row pointer) mapping each cell to its owned face range.
- **Assembly:** `assemble_ldu()` enumerates internal faces analytically. For a cell at `(i,j)`, right-neighbour face has `owner = i*N+j`, `neighbour = i*N+j+1`; top-neighbour face has `owner = i*N+j`, `neighbour = (i+1)*N+j`. Coefficients are `1/h²`. The diagonal accumulates contributions from all connected faces plus boundary ghost-cell treatment.
- **GS sweep:** The `bPrime` push-forward pattern. Initialise `bPrime = b`. For each cell in order: read `psii = bPrime[celli]`, subtract upper-triangle contributions via `psii -= upper[f] * psi[uAddr[f]]` for owned faces, divide by `diag[celli]`, then push lower-triangle contributions forward via `bPrime[uAddr[f]] -= lower[f] * psii`. This avoids a reverse-lookup for lower-triangle neighbours — the same trick OpenFOAM uses.

### DIA Solver (`dia_solver.c`)

Replaces face-indexed indirection with direct strided access:

- **Storage:** Five contiguous arrays of length `N²` — `d0[i]` (main diagonal), `upper1[i]` and `lower1[i]` (offset ±1, horizontal neighbours), `upperN[i]` and `lowerN[i]` (offset ±N, vertical neighbours). No index arrays. No `ownerStart`.
- **Assembly:** `assemble_dia()` fills each diagonal array directly. For cell `idx = i*N+j`: if `j < N-1`, set `upper1[idx] = 1/h²`; if `i < N-1`, set `upperN[idx] = 1/h²`. Boundary cells get zero in the corresponding diagonal entry and adjusted `d0`.
- **GS sweep:** A single cell loop. For cell `i`: `psi[i] = (b[i] - lower1[i]*psi[i-1] - lowerN[i]*psi[i-N] - upper1[i]*psi[i+1] - upperN[i]*psi[i+N]) / d0[i]`. All five array reads are contiguous or fixed-stride — the hardware prefetcher handles them without any indirection.

### Validation

Both solvers must produce cell-identical results. The validation strategy:

1. **Manufactured solution convergence:** After running to convergence, compute L2 error against `u_exact`. Both solvers must report the same error to machine precision.
2. **Cross-solver agreement:** After a fixed number of sweeps (e.g. 100), assert `‖psi_ldu − psi_dia‖∞ < 1e-10`. Any discrepancy indicates a boundary-handling bug.
3. **Residual monitoring:** `compute_residual()` applies the 5-point stencil directly (no solver-specific arrays) — an independent reference both solvers check against.
4. **Convergence order:** Run at N = 16, 32, 64, 128 and verify O(h²) error reduction, confirming the discretisation is correct.

## Project Structure

```
poisson_bench/
├── Makefile
├── README.md
├── include/
│   ├── problem.h          # Grid params, source/exact solution, assembly
│   ├── ldu_solver.h        # LDU struct + gs_sweep_ldu()
│   ├── dia_solver.h        # DIA struct + gs_sweep_dia()
│   └── utils.h             # Residual, L2 error, allocation helpers
├── src/
│   ├── problem.c           # compute_source(), compute_exact(), assemble_ldu(), assemble_dia()
│   ├── ldu_solver.c         # gs_sweep_ldu() — bPrime push-forward pattern
│   ├── dia_solver.c         # gs_sweep_dia() — direct strided access
│   ├── utils.c              # compute_l2_error(), compute_residual()
│   └── main.c              # Parse N, run both solvers, validate, benchmark
└── scripts/
    └── run_likwid.sh        # LIKWID wrapper for all mesh sizes + perf groups
```

## Building

```bash
# Standard build
make

# Build with LIKWID profiling support
make LIKWID=1

# Clean
make clean
```

**Compiler flags:** `-O3 -march=native -ffast-math -Wall`

When built with `LIKWID=1`, adds `-DLIKWID_PERFMON` and links `-llikwid`.

**Requirements:** GCC (tested on GCC 12+), GNU Make. LIKWID for profiling (optional).

## Usage

```bash
# Run both solvers on a 128×128 grid, validate, and print timing
./poisson_bench 128

# Run with a specific sweep count
./poisson_bench 256 --sweeps 500

# LIKWID profiling (requires LIKWID build)
likwid-perfctr -g FLOPS_DP -m ./poisson_bench 512
```

## Profiling

The LIKWID wrapper script runs both solvers across multiple mesh sizes with the performance groups most relevant to the memory-bound hypothesis:

```bash
./scripts/run_likwid.sh
```

This runs groups `FLOPS_DP`, `L2CACHE`, `L3CACHE`, `MEM` at N = 32, 128, 512, 1024. The key metrics to compare between LDU and DIA:

- **Effective memory bandwidth** vs theoretical peak (dual-channel DDR5 ~50–70 GB/s)
- **L3 miss rate** — should drop dramatically for DIA due to prefetchable access patterns
- **FLOP/s utilisation** — should increase from the 0.31% baseline once cores are no longer starved

## Hardware

- **CPU:** AMD Ryzen 7 250 (Zen4), 8C/16T, 3.29 GHz base
- **Memory:** Dual-channel DDR5
- **Profiling:** LIKWID 5.x (groups CACHE, MEM, DATA, L2CACHE, L3CACHE, FLOPS_DP)
- **OS:** Linux (Fedora/Ubuntu)

## What Comes Next

Phase 2 wires the proven DIA kernel into OpenFOAM as a drop-in smoother plugin:

- A new `structuredGaussSeidelSmoother` class inheriting from `lduMatrix::smoother`, registered via `addToRunTimeSelectionTable` — users switch to it by changing one line in `fvSolution`:
  ```
  smoother  structuredGaussSeidel;
  ```
- At solver setup, the constructor detects whether the mesh is structured (uniform `deltaCoeffs`, consistent neighbour counts) and converts the `lduMatrix` into DIA format (one-off cost, amortised over solver iterations). Falls back to default GS if the mesh is unstructured.
- Validation against the default smoother on cavity and other structured tutorials — residual convergence and solution fields must match within solver tolerance.

The kernel is kept clean and portable for a future port to the Uttunga RISC-V posit co-processor board, where the same addressing advantage applies but the arithmetic is posit rather than IEEE 754.

## Author

Amartya — HPC Software Engineer, Calligo Technologies, Bengaluru

## License

Internal / proprietary. Contact Calligo Technologies for licensing enquiries.