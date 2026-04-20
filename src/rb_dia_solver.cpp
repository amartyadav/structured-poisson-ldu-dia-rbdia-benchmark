#include "rb_dia_solver.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void assemble_rbdia(RBDIAMatrix *mat, double *source, int N)
{
    mat->N = N;
    mat->nCells = N * N * N;
    mat->source = (double *)malloc(sizeof(double) * mat->nCells);
    mat->bPrime = (double *)malloc(sizeof(double) * mat->nCells);
    mat->diag = (double *)malloc(sizeof(double) * mat->nCells);
    mat->psi = (double *)calloc(mat->nCells, sizeof(double)); // initial guess for the solver (solution vector) = 0 (hence, calloc);

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
void gs_sweep_rbdia(RBDIAMatrix *mat)
{
    int N = mat->N;

    for (int idx = 0; idx < mat->nCells; idx++)
    {
        mat->bPrime[idx] = mat->source[idx];
    }

    // red sweep
    #pragma omp parallel for
    for(int idx = 0; idx < mat->nCells; idx++)
    {
        // if(omp_get_thread_num() == 0) { printf("Total running threads = %d", omp_get_num_threads());;}

        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        if((i+j+k) % 2 != 0) continue;//skipping half the cells (possible area of improvement/optimisation)

        double psii = mat->bPrime[idx];

        //forward neighbours
        if (j < N - 1) { psii -= -1.0 * mat->psi[idx + 1]; }
        if (i < N - 1) { psii -= -1.0 * mat->psi[idx + N]; }
        if (k < N - 1) { psii -= -1.0 * mat->psi[idx + N * N]; }

        // backward neighbours
        if(j > 0) { psii -= -1.0 * mat->psi[idx - 1]; }
        if(i > 0) { psii -= -1.0 * mat->psi[idx - N]; }
        if(k > 0) { psii -= -1.0 * mat->psi[idx - N * N]; }

        psii /= mat->diag[idx];
        mat->psi[idx] = psii;
    }

    // black sweep
    #pragma omp parallel for
    for(int idx = 0; idx < mat->nCells; idx++)
    {
        // if(omp_get_thread_num() == 0) { printf("Total running threads = %d", omp_get_num_threads());;}

        int k = idx / (N*N);
        int i = (idx / N) % N;
        int j = idx % N;

        if((i + j + k) % 2 != 1) continue;

        double psii = mat->bPrime[idx];
        //forward neighbours
        if (j < N - 1) { psii -= -1.0 * mat->psi[idx + 1]; }
        if (i < N - 1) { psii -= -1.0 * mat->psi[idx + N]; }
        if (k < N - 1) { psii -= -1.0 * mat->psi[idx + N * N]; }

        // backward neighbours
        if(j > 0) { psii -= -1.0 * mat->psi[idx - 1]; }
        if(i > 0) { psii -= -1.0 * mat->psi[idx - N]; }
        if(k > 0) { psii -= -1.0 * mat->psi[idx - N * N]; }

        psii /= mat->diag[idx];
        mat->psi[idx] = psii;
    }
}

// function to enable and store memory tracing for visualisation
void gs_sweep_rbdia_trace(RBDIAMatrix *mat, FILE *fp)
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

        int read_array[3];
        int write_array[3];
        int read_array_idx = 0;
        int write_array_idx = 0;


        double psii = mat->bPrime[idx];

        // upper triangle pull
        // also writing the indices being accessed to the file
        if (j < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + 1];
            read_array[read_array_idx] = idx+1;
            read_array_idx++;
        }
        if (i < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N];
            read_array[read_array_idx] = idx+N;
            read_array_idx++;
        }
        if (k < N - 1)
        {
            psii -= -1.0 * mat->psi[idx + N * N];
            read_array[read_array_idx] = idx + N * N;
            read_array_idx++;
        }

        psii /= mat->diag[idx];

        // lower triangle push
        if (j < N - 1)
        {
            mat->bPrime[idx + 1] -= -1.0 * psii;
            write_array[write_array_idx] = idx+1;
            write_array_idx++;
        }
        if (i < N - 1)
        {
            mat->bPrime[idx + N] -= -1.0 * psii;
            write_array[write_array_idx] = idx+N;
            write_array_idx++;
        }
        if (k < N - 1)
        {
            mat->bPrime[idx + N * N] -= -1.0 * psii;
            write_array[write_array_idx] = idx + N * N;
            write_array_idx++;
        }

        mat->psi[idx] = psii;

        fprintf(fp, "{\"solver\":\"DIA\",\"cell\":%d,\"read\":[", idx);
        for (int r = 0; r < read_array_idx; r++)
        {
            if (r > 0) fprintf(fp, ",");
            fprintf(fp, "%d", read_array[r]);
        }
        fprintf(fp, "],\"write\":[");

        for (int w = 0; w < write_array_idx; w++)
        {
            if (w > 0) fprintf(fp, ",");
            fprintf(fp, "%d", write_array[w]);
        }
        fprintf(fp, "]}\n");
    }
}

void free_rbdia(RBDIAMatrix *mat)
{
    free(mat->bPrime);
    free(mat->diag);
    free(mat->psi);
    free(mat->source);
}
