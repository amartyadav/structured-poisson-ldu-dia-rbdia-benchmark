#include <math.h>
#include <stdio.h>
#include <iostream>
#include "problem.hpp"
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void test_math_sanity();

void test_math_sanity()
{
    // Check 1: M_PI, sin(M_PI/2), sin(M_PI), cos(0)
    std::cout << "=== Math Sanity Check 1: Basic Constants and Functions ===\n";
    std::cout << "M_PI = " << M_PI << " (should be approximately 3.141593)\n";
    std::cout << "sin(M_PI/2) = " << sin(M_PI / 2) << " (should be approximately 1.0)\n";
    std::cout << "sin(M_PI) = " << sin(M_PI)    << " (should be approximately 0.0)\n";
    std::cout << "cos(0) = " << cos(0)          << " (should be approximately 1.0)\n";

    // Check 2: For N=50, compute and print the coordinates of a few known cells — cell 0, cell at the midpoint (nCells/2), and a cell somewhere in the interior. Print i, j, k, x, y, z for each. Check that the midpoint cell gives coordinates near (0.5, 0.5, 0.5). If any coordinate is 0 when it shouldn't be, you have the integer division bug.
    std::cout << "\n=== Math Sanity Check 2: Cell Indexing and Coordinates ===\n";
    int N = 50;
    int nCells = N * N * N;
    int test_indices[] = {0, nCells / 2, nCells / 2 + N * N + N + 1}; // cell 0, midpoint cell, and an interior cell
    for (int idx : test_indices)
    {
        int k = idx / (N * N);
        int i = (idx / N) % N;
        int j = idx % N;
        double h = 1.0 / N;
        double x = (j + 0.5) * h;
        double y = (i + 0.5) * h;
        double z = (k + 0.5) * h;
        std::cout << "Cell index: " << idx << " -> (i, j, k) = (" << i << ", " << j << ", " << k << ") -> (x, y, z) = (" << x << ", " << y << ", " << z << ")\n";
    }

    // Check 3:
    // Print 1 / N, 1.0 / N, (double)1 / N and see which gives 0 vs the correct value.On the Uttunga posit compiler, the cast behaviour might differ from GCC on x86.
    std::cout << "\n=== Math Sanity Check 3: Integer Division Test ===\n";
    std::cout << "1/ N = " << 1 / N << " (should be 0)\n";
    std::cout << "1.0/ N = " << 1.0 / N << " (should be the correct value)\n";
    std::cout << "(double)1/ N = " << (double)1 / N << " (should be the correct value)\n";

    // Check 4: End-to-end source check.
    // Call compute_source for N=4 (only 64 cells) and print every value. Compare against x86 output of the same. The first cell and the last cell where values diverge will tell you exactly which computation is going wrong.
    std::cout << "\n=== Math Sanity Check 4: End-to-End Source Computation ===\n";
    N = 4;
    nCells = N * N * N;
    double *b = (double *)malloc(sizeof(double) * nCells);
    compute_source(b, N);
    for (int idx = 0; idx < nCells; idx++)
    {
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "b[" << idx << "] = " << b[idx] << "\n";
    }
    free(b);
}