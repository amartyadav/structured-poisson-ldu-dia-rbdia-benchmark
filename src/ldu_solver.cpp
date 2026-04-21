// Initiating the setup of the LDU matrix addressing scheme for the solver

#include "ldu_solver.hpp"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/// @brief Assemble LDU storage and addressing for a 3D N x N x N Laplacian grid.
/// @param mat Output matrix container; internal arrays are allocated and initialized.
/// @param source Right-hand-side vector of length N^3.
/// @param N Grid resolution per dimension.
void assemble_ldu(LDUMatrix *mat, double *source, int N)
{
    // Allocate topology and coefficient arrays for the LDU stencil representation.

    mat->nFaces = 3 * N * N * (N - 1); // ((N - 1) * N * N) + (N * (N - 1) * N) + (N * N * (N - 1));
    mat->nCells = N * N * N;
    mat->nFaces = 3 * N * N * (N - 1);
    mat->diag = (double *)calloc(mat->nCells, sizeof(double));
    mat->upper = (double *)malloc(sizeof(double) * mat->nFaces);
    mat->lower = (double *)malloc(sizeof(double) * mat->nFaces);
    mat->uAddr = (int *)malloc(sizeof(int) * mat->nFaces);
    mat->ownerStart = (int *)malloc(sizeof(int) * (mat->nCells + 1));
    mat->psi = (double *)calloc(mat->nCells, sizeof(double)); // initial guess for the solver (solution vector) = 0;
    mat->bPrime = (double *)malloc(sizeof(double) * mat->nCells);
    mat->source = (double *)malloc(sizeof(double) * mat->nCells);

    // Copy source into owned storage so the caller can manage its own memory.
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->source[idx] = source[idx];
    }

    // Running face counter while constructing owner->upper addressing.
    int f = 0;
    mat->ownerStart[0] = 0;

    for (int idx = 0; idx < N * N * N; idx++)
    {
        // Flattened cell index -> (k, i, j) coordinates.
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        // Start offset of faces owned by this cell.
        mat->ownerStart[idx] = f;

        // Own +j face if a right neighbor exists.
        if (j < N - 1)
        {
            // Face couples current owner cell (idx) with upper cell (idx + 1).
            mat->lower[f] = -1.0;
            mat->upper[f] = -1.0;
            mat->uAddr[f] = idx + 1;
            mat->diag[idx] += 1.0;
            mat->diag[idx + 1] += 1.0;
            f++;
        }

        // Own +i face if a top neighbor exists.
        if (i < N - 1)
        {
            // Face couples current owner cell (idx) with upper cell (idx + N).
            mat->lower[f] = -1.0;
            mat->upper[f] = -1.0;
            mat->uAddr[f] = idx + N;
            mat->diag[idx] += 1.0;
            mat->diag[idx + N] += 1.0;
            f++;
        }

        // Own +k face if a front/back neighbor exists.
        if (k < N - 1)
        {
            mat->lower[f] = -1.0;
            mat->upper[f] = -1.0;
            mat->uAddr[f] = idx + (N * N);
            mat->diag[idx] += 1.0;
            mat->diag[idx + (N * N)] += 1.0;
            f++;
        }
    }

    // Dirichlet boundary contribution: add boundary-face weights to the diagonal.
    for (int idx = 0; idx < N * N * N; idx++)
    {
        // Flattened cell index -> (k, i, j) coordinates.
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        if (j == 0)
        {
            mat->diag[idx] += 2.0;
        }
        if (j == N - 1)
        {
            mat->diag[idx] += 2.0;
        }
        if (i == 0)
        {
            mat->diag[idx] += 2.0;
        }
        if (i == N - 1)
        {
            mat->diag[idx] += 2.0;
        }
        if (k == 0)
        {
            mat->diag[idx] += 2.0;
        }
        if (k == N - 1)
        {
            mat->diag[idx] += 2.0;
        }
    }

    // Sentinel: end of owned-face range for the last cell.
    mat->ownerStart[N * N * N] = f; // should this be mat->nCells = f  instead?
    printf("f = %d, nFaces = %d\n", f, mat->nFaces);
    // Sanity check that constructed addressing matches expected face count.
    assert(f == mat->nFaces);
}

/// @brief Perform one forward Gauss-Seidel sweep using LDU owner/upper addressing.
/// @param mat Matrix and vectors for the current nonlinear/linear iteration state.
void gs_sweep_ldu(LDUMatrix *mat)
{
    // Reset transformed RHS at the beginning of each sweep.
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        // Face range owned by cell idx is [f_start, f_end).
        int f_start = mat->ownerStart[idx];
        int f_end = mat->ownerStart[idx + 1];

        // Start from current RHS contribution for this row.
        double psii = mat->bPrime[idx];

        // Upper-triangular pull: subtract stale values from not-yet-updated neighbors.
        for (int f = f_start; f < f_end; f++)
        {
            psii -= mat->upper[f] * mat->psi[mat->uAddr[f]];
        }

        // Apply diagonal scaling once per row.
        psii /= mat->diag[idx];

        // Lower-triangular push: inject updated psi_i into future RHS entries.
        for (int f = f_start; f < f_end; f++)
        {
            mat->bPrime[mat->uAddr[f]] -= mat->lower[f] * psii;
        }

        // Commit updated solution component.
        mat->psi[idx] = psii;
    }
}

/// @brief Release all heap-allocated arrays owned by the LDU matrix container.
/// @param mat Matrix whose internal buffers are deallocated.
void free_ldu(LDUMatrix *mat)
{
    free(mat->bPrime);
    free(mat->diag);
    free(mat->lower);
    free(mat->upper);
    free(mat->ownerStart);
    free(mat->psi);
    free(mat->source);
    free(mat->uAddr);
}
