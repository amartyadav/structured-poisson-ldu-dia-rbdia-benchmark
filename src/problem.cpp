#include "problem.hpp"
#include <math.h>

#define M_PI 3.14159265358979323846

/// @brief Build RHS for the manufactured Poisson problem on an N x N x N grid.
/// @param b Output source vector of length N^3.
/// @param N Grid resolution per dimension.
void compute_source(double *b, int N)
{
    // Uniform cell-centered spacing in the unit cube.
    double h = 1.0 / N;
    double h_square = h * h;
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Flattened index and cell-center coordinates.
                int idx = k * N * N + i * N + j;
                double x = (j + 0.5) * h;
                double y = (i + 0.5) * h;
                double z = (k + 0.5) * h;

                // Discrete RHS for u = sin(pi x) sin(pi y) sin(pi z), scaled by h^2.
                b[idx] = h_square * (3.0 * M_PI * M_PI) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
            }
        }
    }
}

/// @brief Compute manufactured exact solution at cell centers.
/// @param u_exact Output exact solution vector of length N^3.
/// @param N Grid resolution per dimension.
void compute_exact(double *u_exact, int N)
{
    // Uniform cell-centered spacing in the unit cube.
    double h = 1.0 / N;
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Flattened index and cell-center coordinates.
                int idx = k * N * N + i * N + j;
                double x = (j + 0.5) * h;
                double y = (i + 0.5) * h;
                double z = (k + 0.5) * h;

                // Manufactured exact field.
                u_exact[idx] = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
            }
        }
    }
}

/// @brief Compute RMS L2 error between numerical and exact solutions.
/// @param u Numerical solution vector.
/// @param u_exact Exact reference solution vector.
/// @param N Grid resolution per dimension.
/// @return Root-mean-square error over all N^3 cells.
double compute_l2_error(double *u, double *u_exact, int N)
{
    double sum = 0.0;
    int nCells = N * N * N;
    for (int i = 0; i < nCells; i++)
    {
        // Accumulate squared pointwise error.
        sum += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
    }
    return sqrt(sum / nCells);
}

/// @brief Compute RMS residual norm ||b - A u|| on the 3D structured stencil.
/// @param u Current solution iterate.
/// @param b Right-hand-side vector.
/// @param N Grid resolution per dimension.
/// @return Root-mean-square residual over all N^3 cells.
double compute_residual(double *u, double *b, int N)
{
    double sum = 0.0;
    int nCells = N * N * N;
    for (int idx = 0; idx < nCells; idx++)
    {
        // Flattened index -> (k, i, j) coordinates.
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;

        // Interior diagonal plus boundary-face Dirichlet contributions.
        double diag = 6.0;
        if (j == 0)
        {
            diag += 1.0;
        }
        if (j == N - 1)
        {
            diag += 1.0;
        }
        if (i == 0)
        {
            diag += 1.0;
        }
        if (i == N - 1)
        {
            diag += 1.0;
        }
        if (k == 0)
        {
            diag += 1.0;
        }
        if (k == N - 1)
        {
            diag += 1.0;
        }

        // Apply local stencil row A(idx, : ) to u.
        double Au_idx = diag * u[idx];
        if (j < N - 1)
            Au_idx -= u[idx + 1];
        if (j > 0)
            Au_idx -= u[idx - 1];
        if (i < N - 1)
            Au_idx -= u[idx + N];
        if (i > 0)
            Au_idx -= u[idx - N];
        if (k < N - 1)
            Au_idx -= u[idx + N * N];
        if (k > 0)
            Au_idx -= u[idx - N * N];

        // Accumulate squared residual.
        sum += (b[idx] - Au_idx) * (b[idx] - Au_idx);
    }
    return sqrt(sum / nCells);
}
