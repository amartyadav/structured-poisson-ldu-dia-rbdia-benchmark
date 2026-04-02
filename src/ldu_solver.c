// Initiating the setup of the LDU matrix addressing scheme for the solver

#include "ldu_solver.h"
#include <stdlib.h>

void assemble_ldu(LDUMatrix *mat, double *source, int N)
{
    // allocating the arrays

    mat->nCells = N * N;
    mat->nFaces = 2 * N * (N-1);
    mat->diag = calloc(mat->nCells, sizeof(double));
    mat->upper = malloc(sizeof(double) * mat->nFaces);
    mat->lower = malloc(sizeof(double) * mat->nFaces);
    mat->uAddr = malloc(sizeof(int) * mat->nFaces);
    mat->ownerStart = malloc(sizeof(int) * (mat->nCells + 1));
    mat->psi = calloc(mat->nCells, sizeof(double)); // initial guess for the solver (solution vector) = 0;
    mat->bPrime = malloc(sizeof(double) * mat->nCells);
    mat->source = malloc(sizeof(double) * mat->nCells);

    // copying the source
    for(int idx = 0; idx < mat->nCells; idx++)
    {
        mat->source[idx] = source[idx];
    }

    int f = 0;
    mat->ownerStart[0] = 0;

    for(int idx = 0; idx < N * N; idx++)
    {
        int i = idx / N; // row
        int j = idx % N; // column

        // this cells face starts from here
        mat->ownerStart[idx] = f;

        // does this cell (current cell) own a right face (have a right neighbour) ?
        if (j < N - 1)
        {
            // there is a face between idx (cell) and idx + 1 (cell + 1)
            mat->lower[f] = -1.0;
            mat->upper[f] = -1.0;
            mat->uAddr[f] = idx+1;
            mat->diag[idx] += 1.0;
            mat->diag[idx+1] += 1.0;
            f++;
        }

        // does this cell (current cell) own a top face (have a top neighbour) ?
        if (i < N - 1)
        {
            // there is a cell between idx (cell) and idx + N (cell + N)
            mat->lower[f] = -1.0;
            mat->upper[f] = -1.0;
            mat->uAddr[f] = idx+N;
            mat->diag[idx] += 1.0;
            mat->diag[idx+N] += 1.0;
            f++;
        }
    }

    // After the face-building loop, fix up boundary cells for Dirichlet BC
    for (int idx = 0; idx < N*N; idx++) {
        int i = idx / N;
        int j = idx % N;
        if (j == 0)     { mat->diag[idx] += 1.0; /* source += 0 since sin(0)=0 */ }
        if (j == N-1)   { mat->diag[idx] += 1.0; }
        if (i == 0)     { mat->diag[idx] += 1.0; }
        if (i == N-1)   { mat->diag[idx] += 1.0; }
    }

    mat->ownerStart[N*N] = f; //should this be mat->nCells = f  instead?
}

void gs_sweep_ldu(LDUMatrix *mat)
{
    // reset bPrime to source at the beginning of every sweep
    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    for(int idx = 0; idx < mat->nCells; idx++)
    {
        int f_start = mat->ownerStart[idx];
        int f_end = mat->ownerStart[idx+1];

        double psii = mat->bPrime[idx];

        // upper triangle pull -> subtracting the stale neighbours from the unsolved (upper triangle)
        for(int f = f_start; f < f_end; f++)
        {
            psii -= mat->upper[f] * mat->psi[mat->uAddr[f]];
        }

        // dividing by the diagonal (pivot) -> strictly needs to be outside the loop (otherwise div. can happen twice or thrice or once depending on the number of faces each cell owns)
        psii /=mat->diag[idx];

        // lower triangle push -> pushing the updated values into future cells' RHS
        for(int f = f_start; f < f_end; f++)
        {
            mat->bPrime[mat->uAddr[f]] -= mat->lower[f] * psii;
        }

        // committing the updated value to psii
        mat->psi[idx] = psii;
    }
}



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
