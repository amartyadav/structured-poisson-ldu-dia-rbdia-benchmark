#include "problem.h"
#include <math.h>

#define M_PI 3.14159265358979323846

void compute_source(double *b, int N)
{
    int idx = 0;
    double x, y;
    double h = 1.0/N;
    double h_square = h * h;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            idx = i * N + j;
            x = (j + 0.5) * h; // x = (j + 0.5) * h (h = 1 / N)
            y = (i + 0.5) * h; // y = (i + 0.5) * h (h = 1 / N)
            b[idx] = h_square * (2*M_PI*M_PI) * sin(M_PI*x) * sin(M_PI*y);
        }
    }
}

void compute_exact(double *u_exact, int N)
{
    int idx = 0;
    double x, y;
    double h = 1.0 / N;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            idx = i * N + j;
            x = (j + 0.5) * h; // same as above
            y = (i + 0.5) * h; // same as above
            u_exact[idx] = sin(M_PI*x) * sin(M_PI*y);
        }
    }
}

double compute_l2_error(double *u, double *u_exact, int N)
{
    double sum = 0.0;
    for(int i = 0; i < N*N; i++) {
        sum += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
    }
    return sqrt(sum / (N * N));
}

// Computes the RMS residual ‖b − A·u‖ / N for the discrete 2D Poisson
// system with homogeneous Dirichlet BCs (diagonal = 4 for all cells).
// For cell (i,j) with idx = i*N+j:
//   Au_idx = 4·u[idx]         ← diagonal (4 neighbours or Dirichlet ghost)
//          − u[idx±1]         ← left/right if interior
//          − u[idx±N]         ← bottom/top if interior
//   r[idx] = b[idx] − Au_idx
//   return sqrt( Σ r² / N² )
double compute_residual(double *u, double *b, int N)
{
    double sum = 0.0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            int idx = i*N+j;
            double Au_idx = 4.0 * u[idx];  // diagonal is always 4
            if (j < N - 1) Au_idx -= u[idx+1]; else Au_idx -= 0.0; // Dirichlet = 0
            if (j > 0)     Au_idx -= u[idx-1]; else Au_idx -= 0.0;
            if (i < N - 1) Au_idx -= u[idx+N]; else Au_idx -= 0.0;
            if (i > 0)     Au_idx -= u[idx-N]; else Au_idx -= 0.0;
            sum += (b[idx] - Au_idx) * (b[idx] - Au_idx);
        }
    }
    return sqrt(sum / (N*N));
}
