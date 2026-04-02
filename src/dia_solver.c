#include "dia_solver.h"
#include <stdlib.h>

void assemble_dia(DIAMatrix *mat, double *source, int N)
{
    mat->N = N;
    mat->nCells = N * N;
    mat->source = malloc(sizeof(double) * mat->nCells);
    mat->bPrime = malloc(sizeof(double) * mat->nCells);
    mat->diag = malloc(sizeof(double) * mat->nCells);
    mat->psi = calloc(mat->nCells, sizeof(double)); // initial guess for the solver (solution vector) = 0 (hence, calloc);

    // copying the source
    for(int idx = 0; idx < mat->nCells; idx++)
    {
        mat->source[idx] = source[idx];
    }

    // diagonal coefficients: 4.0 for every cell (Dirichlet BCs already baked in)
    for(int cell = 0; cell < mat->nCells; cell++)
    {
        mat->diag[cell] = 4.0;
    }
}

void gs_sweep_dia(DIAMatrix *mat)
{
    // reset bPrime to source at the beginning of every sweep
    for(int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    // interior sweep : i = 1 -> N - 2, j = 1 -> N - 2;
    for(int i = 0; i < mat->N - 1; i++)
    {
        // interior cells (all neighbours (top and right) exist)
        for(int j = 0; j < mat-> N - 1; j++)
        {
            int idx =i * mat->N + j;
            double psii = mat->bPrime[idx];

            // upper triangle pull -> subtracting the stale values from the unsolved (upper) neighbours
            psii -= -1.0 * mat->psi[idx + 1];
            psii -= -1.0 * mat->psi[idx+mat->N];

            // dividing by the diagnol (pivot)
            psii /= mat->diag[idx];

            // lower triangle push -> pushing the updated values into the forward neighbours' bPrime
            mat->bPrime[idx + 1] -= -1.0 * psii;
            mat->bPrime[idx + mat->N] -= -1.0 * psii;

            // commit the updated value into psi
            mat->psi[idx] = psii;
        }

        // last column of this row i (only top neighbour exists)
        int idx = i * mat->N + (mat->N - 1);
        double psii = mat->bPrime[idx];
        psii -= -1.0 * mat->psi[idx + mat->N];
        psii /= mat->diag[idx];

        mat->bPrime[idx + mat->N] -= -1.0 * psii;
        mat->psi[idx] = psii;
    }

    // last row (i = N - 1); no top neighbour (only right neighbour exists)
    for(int j = 0; j < mat->N - 1; j++)
    {
    int idx = (mat->N - 1) * mat->N + j;
    double psii = mat->bPrime[idx];
    psii -= -1.0 * mat->psi[idx + 1];
    psii /= mat->diag[idx];

    mat->bPrime[idx + 1] -= -1.0 * psii;
    mat->psi[idx] = psii;
    }

    // very last cell (no neighbours in either direction)
    int idx = mat->nCells - 1;
    double psii = mat->bPrime[idx];
    psii /= mat->diag[idx];
    mat->psi[idx] = psii;
}

void free_dia(DIAMatrix *mat)
{
    free(mat->bPrime);
    free(mat->diag);
    free(mat->psi);
    free(mat->source);
}
