#include "dia_solver.hpp"
#include <stdio.h>
#include <stdlib.h>

/// @brief Assemble DIA storage for a 3D Laplacian system on an N x N x N grid.
/// @param mat Output matrix container; internal arrays are allocated and initialized.
/// @param source Right-hand-side vector of length N^3.
/// @param N Grid resolution per dimension.
void assemble_dia(DIAMatrix *mat, double *source, int N)
{
    // Store mesh metadata used by all sweep routines.
    mat->N = N;
    mat->nCells = N * N * N;

    // Allocate contiguous arrays for source, transformed RHS, diagonal, and solution.
    mat->source = (double *)malloc(sizeof(double) * mat->nCells);
    mat->bPrime = (double *)malloc(sizeof(double) * mat->nCells);
    mat->diag = (double *)malloc(sizeof(double) * mat->nCells);
    mat->psi = (double *)calloc(mat->nCells, sizeof(double)); // Initial guess is zero.

    // Copy source into owned storage so caller memory can be independent.
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->source[idx] = source[idx];
    }

    // // diagonal coefficients: 6.0 for every cell (Dirichlet BCs already baked in) (changed from 4.0 as we have 6 neighbours instead of 4 for the diagonol cells)
    // for(int cell = 0; cell < mat->nCells; cell++)
    // {
    //     mat->diag[cell] = 6.0;
    // }

    for (int cell = 0; cell < mat->nCells; cell++)
    {
        // Flattened index -> (k, i, j) coordinates.
        int k = cell / (N * N);
        int i = (cell / N) % N;
        int j = cell % N;

        // Interior 7-point stencil diagonal is 6.0.
        // Boundary planes increase the diagonal to encode Dirichlet treatment.
        mat->diag[cell] = 6.0;
        if (j == 0)
        {
            mat->diag[cell] += 1.0;
        }
        if (j == N - 1)
        {
            mat->diag[cell] += 1.0;
        }
        if (i == 0)
        {
            mat->diag[cell] += 1.0;
        }
        if (i == N - 1)
        {
            mat->diag[cell] += 1.0;
        }
        if (k == 0)
        {
            mat->diag[cell] += 1.0;
        }
        if (k == N - 1)
        {
            mat->diag[cell] += 1.0;
        }
    }
}

/// @brief Perform one forward Gauss-Seidel sweep using DIA-style neighbor access.
/// @param mat Matrix and vectors for the current nonlinear/linear iteration state.
void gs_sweep_dia(DIAMatrix *mat)
{
    int N = mat->N;

    // Reset transformed RHS from the original source each sweep.
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        // Flattened index -> (k, i, j) coordinates.
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        // Start from current RHS contribution for this row.
        double psii = mat->bPrime[idx];

        // Upper-triangular pull: read not-yet-updated forward neighbors.
        if (j < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + 1];
        }
        if (i < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N];
        }
        if (k < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N * N];
        }

        // Diagonal scaling for the new solution value.
        psii /= mat->diag[idx];

        // Lower-triangular push: propagate effect of psi_i into later bPrime entries.
        if (j < N - 1)
        {
            mat->bPrime[idx + 1] -= -1.0 * psii;
        }
        if (i < N - 1)
        {
            mat->bPrime[idx + N] -= -1.0 * psii;
        }
        if (k < N - 1)
        {
            mat->bPrime[idx + N * N] -= -1.0 * psii;
        }

        // Store updated solution component.
        mat->psi[idx] = psii;
    }
}

/// @brief Gauss-Seidel sweep with read/write index tracing for visualization/profiling.
/// @param mat Matrix and vectors for the current iteration state.
/// @param fp Output stream receiving one JSON line per processed cell.
void gs_sweep_dia_trace(DIAMatrix *mat, FILE *fp)
{
    int N = mat->N;

    // Reset transformed RHS from source before this sweep.
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        // Flattened index -> (k, i, j) coordinates.
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        // At most three neighbors in +j, +i, +k directions are touched.
        int read_array[3];
        int write_array[3];
        int read_array_idx = 0;
        int write_array_idx = 0;

        double psii = mat->bPrime[idx];

        // Upper-triangular pull and read-index capture.
        if (j < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + 1];
            read_array[read_array_idx] = idx + 1;
            read_array_idx++;
        }
        if (i < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N];
            read_array[read_array_idx] = idx + N;
            read_array_idx++;
        }
        if (k < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N * N];
            read_array[read_array_idx] = idx + N * N;
            read_array_idx++;
        }

        psii /= mat->diag[idx];

        // Lower-triangular push and write-index capture.
        if (j < N - 1)
        {
            mat->bPrime[idx + 1] -= -1.0 * psii;
            write_array[write_array_idx] = idx + 1;
            write_array_idx++;
        }
        if (i < N - 1)
        {
            mat->bPrime[idx + N] -= -1.0 * psii;
            write_array[write_array_idx] = idx + N;
            write_array_idx++;
        }
        if (k < N - 1)
        {
            mat->bPrime[idx + N * N] -= -1.0 * psii;
            write_array[write_array_idx] = idx + N * N;
            write_array_idx++;
        }

        mat->psi[idx] = psii;

        // Emit JSON trace for this cell: solver tag, cell index, read list, write list.
        fprintf(fp, "{\"solver\":\"DIA\",\"cell\":%d,\"read\":[", idx);
        for (int r = 0; r < read_array_idx; r++)
        {
            if (r > 0)
                fprintf(fp, ",");
            fprintf(fp, "%d", read_array[r]);
        }
        fprintf(fp, "],\"write\":[");

        for (int w = 0; w < write_array_idx; w++)
        {
            if (w > 0)
                fprintf(fp, ",");
            fprintf(fp, "%d", write_array[w]);
        }
        fprintf(fp, "]}\n");
    }
}

/// @brief Release all heap-allocated arrays owned by the DIA matrix container.
/// @param mat Matrix whose internal buffers are deallocated.
void free_dia(DIAMatrix *mat)
{
    free(mat->bPrime);
    free(mat->diag);
    free(mat->psi);
    free(mat->source);
}
