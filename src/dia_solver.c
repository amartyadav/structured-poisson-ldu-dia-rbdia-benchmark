#include "dia_solver.h"
#include <stdlib.h>

void assemble_dia(DIAMatrix *mat, double *source, int N)
{
    mat->N = N;
    mat->nCells = N * N * N;
    mat->source = malloc(sizeof(double) * mat->nCells);
    mat->bPrime = malloc(sizeof(double) * mat->nCells);
    mat->diag = malloc(sizeof(double) * mat->nCells);
    mat->psi = calloc(mat->nCells, sizeof(double)); // initial guess for the solver (solution vector) = 0 (hence, calloc);

    // copying the source
    for(int idx = 0; idx < mat->nCells; idx++)
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
        int k = cell / (N * N);
        int i = (cell / N) % N;
        int j = cell % N;
    
        mat->diag[cell] = 6.0;
        if (j == 0)     { mat->diag[cell] += 1.0; }
        if (j == N - 1) { mat->diag[cell] += 1.0; }
        if (i == 0)     { mat->diag[cell] += 1.0; }
        if (i == N - 1) { mat->diag[cell] += 1.0; }
        if (k == 0)     { mat->diag[cell] += 1.0; }
        if (k == N - 1) { mat->diag[cell] += 1.0; }
    }
}

// updated gs_sweep_dia for 3d version, without the explicit handling for edge cases.
void gs_sweep_dia(DIAMatrix *mat)
{
    int N = mat->N;

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        double psii = mat->bPrime[idx];

        // upper triangle pull
        if (j < N - 1) { psii -= -1.0 * mat->psi[idx + 1]; }
        if (i < N - 1) { psii -= -1.0 * mat->psi[idx + N]; }
        if (k < N - 1) { psii -= -1.0 * mat->psi[idx + N * N]; }

        psii /= mat->diag[idx];

        // lower triangle push
        if (j < N - 1) { mat->bPrime[idx + 1]     -= -1.0 * psii; }
        if (i < N - 1) { mat->bPrime[idx + N]      -= -1.0 * psii; }
        if (k < N - 1) { mat->bPrime[idx + N * N]  -= -1.0 * psii; }

        mat->psi[idx] = psii;
    }
}

void free_dia(DIAMatrix *mat)
{
    free(mat->bPrime);
    free(mat->diag);
    free(mat->psi);
    free(mat->source);
}
