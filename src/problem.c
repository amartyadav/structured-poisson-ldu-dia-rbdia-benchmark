#include "problem.h"
#include <math.h>

void compute_source(double *b, int N)
{
    double h = 1.0 / N;
    double h_square = h * h;
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                int idx = k * N * N + i * N + j;
                double x = (j + 0.5) * h;
                double y = (i + 0.5) * h;
                double z = (k + 0.5) * h;
                b[idx] = h_square * (3.0 * M_PI * M_PI) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
            }
        }
    }
}

void compute_exact(double *u_exact, int N)
{
    double h = 1.0 / N;
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                int idx = k * N * N + i * N + j;
                double x = (j + 0.5) * h;
                double y = (i + 0.5) * h;
                double z = (k + 0.5) * h;
                u_exact[idx] = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
            }
        }
    }
}

double compute_l2_error(double *u, double *u_exact, int N)
{
    double sum = 0.0;
    int nCells = N * N * N;
    for (int i = 0; i < nCells; i++)
    {
        sum += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
    }
    return sqrt(sum / nCells);
}

double compute_residual(double *u, double *b, int N)
{
    double sum = 0.0;
    int nCells = N * N * N;
    for (int idx = 0; idx < nCells; idx++)
    {
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        double Au_idx = 6.0 * u[idx];

        if (j < N - 1) Au_idx -= u[idx + 1];
        if (j > 0)     Au_idx -= u[idx - 1];
        if (i < N - 1) Au_idx -= u[idx + N];
        if (i > 0)     Au_idx -= u[idx - N];
        if (k < N - 1) Au_idx -= u[idx + N * N];
        if (k > 0)     Au_idx -= u[idx - N * N];

        sum += (b[idx] - Au_idx) * (b[idx] - Au_idx);
    }
    return sqrt(sum / nCells);
}