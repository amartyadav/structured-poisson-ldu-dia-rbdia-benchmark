#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "problem.hpp"
#include "ldu_solver.hpp"
#include "dia_solver.hpp"
#include "rb_dia_solver.hpp"
#include <chrono>
#include "math_sanity_check.hpp"

#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(a)
#define LIKWID_MARKER_STOP(a)
#define LIKWID_MARKER_CLOSE
#endif

/// @brief Optional legacy multi-size test hook (currently disabled in this file).
void test_multiple_N(void);
/// @brief Run fixed-count sweep profiling for LDU, DIA, and RBDIA solvers.
/// @param N Grid resolution per dimension.
/// @param num_sweeps Number of GS sweeps to execute per solver.
void profile_sweeps(int N, int num_sweeps);
/// @brief Run DIA-only convergence test to a residual tolerance.
/// @param N Grid resolution per dimension.
/// @param tol Residual tolerance.
/// @param max_sweeps Maximum sweep budget.
void convergence_test_LDU(int N, double tol, int max_sweeps);
void convergence_test_DIA(int N, double tol, int max_sweeps);
/// @brief Run RBDIA-only convergence test to a residual tolerance. 
/// @param N Grid resolution per dimension.
/// @param tol Residual tolerance.
/// @param max_sweeps Maximum sweep budget.
void convergence_test_RBDIA(int N, double tol, int max_sweeps);
/// @brief Run DIA and RBDIA convergence tests back-to-back.
/// @param N Grid resolution per dimension.
/// @param tol Residual tolerance.
/// @param max_sweeps Maximum sweep budget.
void convergence_test_combined(int N, double tol, int max_sweeps);

// Helper functions to read input with defaults
static int read_int(const char *prompt, int default_val)
{
    std::string line;
    printf("%s [default=%d]: ", prompt, default_val);
    std::getline(std::cin, line);
    if (line.empty())
        return default_val;
    return std::stoi(line);
}

static double read_double(const char *prompt, double default_val)
{
    std::string line;
    printf("%s [default=%.1e]: ", prompt, default_val);
    std::getline(std::cin, line);
    if (line.empty())
        return default_val;
    return std::stod(line);
}

/// @brief CLI entry point selecting profiling or convergence workflows.
/// @return Exit code (0 on normal completion).
int main()
{
    std::cout << "=== Select your case === \n";
    std::cout << "1. Profile mode: Run a profiling session for a given N and number of sweeps. 1st arg = N, 2nd arg = num_sweeps.\n";
    std::cout << "2. Convergence mode: Run a convergence test for both DIA and RBDIA for a given N, tolerance, and max sweeps. Set OMP_NUM_THREADS=n for the OpenMP threads for RB DIA\n";
    std::cout << "3. DIA convergence mode: Run a convergence test for DIA only for a given N, tolerance, and max sweeps.\n";
    std::cout << "4. RBDIA convergence mode: Run a convergence test for RBDIA only for a given N, tolerance, and max sweeps. Set OMP_NUM_THREADS=n for the OpenMP threads for RB DIA\n";
    std::cout << "5. LDU convergence mode: Run a convergence test for LDU only for a given N, tolerance, and max sweeps.\n";
    std::cout << "6. Math sanity check: Check the sanity and correctness of various math functions, calculations, etc. relevant to this benchmark for x86 and RISC-V Uttunga.\n";
    std::cout << "7. Help: Print usage information.\n";

    std::cout << "Enter the case number (1-7): ";
    int case_number;
    std::cin >> case_number;
    std::cin.ignore(); // consume leftover newline so getline works

    if (case_number == 1)
    {
        int N = read_int("Enter N (grid resolution per dimension)", 50);
        int num_sweeps = read_int("Enter number of sweeps", 1000);
        profile_sweeps(N, num_sweeps);
    }
    else if (case_number == 2)
    {
        int N = read_int("Enter N (grid resolution per dimension)", 50);
        double tol = read_double("Enter residual tolerance", 1e-16);
        int max_sweeps = read_int("Enter maximum number of sweeps", 50000);
        convergence_test_combined(N, tol, max_sweeps);
    }
    else if (case_number == 3)
    {
        int N = read_int("Enter N (grid resolution per dimension)", 50);
        double tol = read_double("Enter residual tolerance", 1e-16);
        int max_sweeps = read_int("Enter maximum number of sweeps", 50000);
        convergence_test_DIA(N, tol, max_sweeps);
    }
    else if (case_number == 4)
    {
        int N = read_int("Enter N (grid resolution per dimension)", 50);
        double tol = read_double("Enter residual tolerance", 1e-16);
        int max_sweeps = read_int("Enter maximum number of sweeps", 50000);
        convergence_test_RBDIA(N, tol, max_sweeps);
    }
    else if (case_number == 5)
    {
        int N = read_int("Enter N (grid resolution per dimension)", 50);
        double tol = read_double("Enter residual tolerance", 1e-16);
        int max_sweeps = read_int("Enter maximum number of sweeps", 50000);
        convergence_test_LDU(N, tol, max_sweeps);
    }
    else if (case_number == 6)
    {
        std::cout << "Running math sanity checks...\n";
        test_math_sanity();
        std::cout << "Math sanity checks completed.\n";
    }
    else if (case_number == 7)
    {
        // Help/usage output when no valid mode is provided.
        printf("\n\n=== USAGE ===\n");
        printf("./poisson_benchmark [mode] [N] [tolerance] [max_sweeps]\n");

        printf("--- mode           - The mode of the benchmark/run that you require.\n");
        printf("\n=== MODES ===\n");
        printf("--- profile             -     !!! ONLY RUN THIS MODE IF COMPILED WITH 'make profile' - Read 'Profiling' section below for details!!!\n");
        printf("                              Run a profiling session for a given N and number of sweeps. 1st arg = N, 2nd arg = num_sweeps.\n");
        printf("--- convergence         -     Run a convergence test for both DIA and RBDIA for a given N, tolerance, and max sweeps. Set OMP_NUM_THREADS=n for the OpenMP threads for RB DIA\n");
        printf("--- diaconvergence      -     Run a convergence test for DIA only for a given N, tolerance, and max sweeps.\n");
        printf("--- rbdiaconvergence    -     Run a convergence test for RBDIA only for a given N, tolerance, and max sweeps. Set OMP_NUM_THREADS=n for the OpenMP threads for RB DIA\n");
        printf("--- lduconvergence      -     Run a convergence test for LDU only for a given N, tolerance, and max sweeps.\n");
        printf("--- mathsanity          -     Check the sanity and correctness of various math functions, calculations, etc. relevant to this benchmark for x86 and RISC-V Uttunga.\n");

        printf("\n=== DEFAULTS ===\n");
        printf("N = 50, tolerance = 1e-16, max_sweeps = 50000\n");
        printf("Press Enter at any prompt to accept the default value.\n");

        printf("\n=== N ===\n");
        printf("The number of cells in each direction. Total cells = N^3. Default is 50.\n");

        printf("\n=== tolerance ===\n");
        printf("The residual tolerance for convergence tests. Default is 1e-16.\n");

        printf("\n=== max_sweeps ===\n");
        printf("The maximum number of sweeps to perform in convergence tests before giving up. Default is 50000.\n");

        printf("\n=== PROFILING ===\n");
        printf("'profile' mode will run a profiling session for a given N and number of sweeps. It will print the working set size for both DIA and LDU, and then run the specified number of sweeps for each solver while timing them with LIKWID markers. The residual after the sweeps will also be printed.\n\n");
        printf("NOTE: 'profile' mode should only be run if the code is compiled with 'make profile', which enables LIKWID performance monitoring. Running 'profile' mode without LIKWID support will not give meaningful results. Also, profiling mode creates three binaries - O1, O2, and O3. Run accordingly.\n");
    }
    else
    {
        std::cerr << "Invalid case number. Please run the program again and select a number between 1 and 6.\n";
        return 1;
    }
    return 0;
}

/// @brief Execute a profiling campaign for all solver variants at fixed sweep count.
/// @param N Grid resolution per dimension.
/// @param num_sweeps Number of sweeps to run for each solver.
void profile_sweeps(int N, int num_sweeps)
{
    // Problem size and topology summary for memory estimates.
    int nCells = N * N * N;
    int nFaces = 3 * N * N * (N - 1);

    printf("=== Profiling run: N=%d, sweeps=%d ===\n", N, num_sweeps);
    printf("Working set (DIA): ~%.1f MB\n",
           (4.0 * nCells * sizeof(double)) / (1024.0 * 1024.0));
    printf("Working set (LDU): ~%.1f MB\n",
           (4.0 * nCells * sizeof(double) + 3.0 * nFaces * sizeof(double) + nFaces * sizeof(int) + (nCells + 1) * sizeof(int)) / (1024.0 * 1024.0));

    // Common manufactured RHS used by all solver variants.
    double *b = (double *)malloc(sizeof(double) * nCells);
    compute_source(b, N);

    // LDU sweep loop with marker region and per-call wall-clock accounting.
    LDUMatrix ldu;
    assemble_ldu(&ldu, b, N);
    printf("\n--- LDU solver: %d sweeps ---\n", num_sweeps);

    LIKWID_MARKER_START("LDU_sweep");
    for (int s = 0; s < num_sweeps; s++)
    {
        static double totalSweepTimeLDU = 0.0;
        static int totalCallsLDU = 0;
        auto t0LDU = std::chrono::high_resolution_clock::now();
        gs_sweep_ldu(&ldu);
        auto t1LDU = std::chrono::high_resolution_clock::now();
        totalSweepTimeLDU += std::chrono::duration<double>(t1LDU - t0LDU).count();
        totalCallsLDU++;
        if (totalCallsLDU % 1000 == 0)
        {
            std::cout << "LDU_sweep: totalTime=" << totalSweepTimeLDU
                      << "s calls=" << totalCallsLDU << std::endl;
        }
    }
    LIKWID_MARKER_STOP("LDU_sweep");

    double res_ldu = compute_residual(ldu.psi, b, N);
    printf("LDU residual after %d sweeps: %.6e\n", num_sweeps, res_ldu);
    free_ldu(&ldu);

    // DIA sweep loop with marker region and aggregate timing.
    DIAMatrix dia;
    assemble_dia(&dia, b, N);
    printf("\n--- DIA solver: %d sweeps ---\n", num_sweeps);

    static double totalSweepTimeDIA = 0.0;
    static int totalCallsDIA = 0;
    auto t0DIA = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("DIA_sweep");
    for (int s = 0; s < num_sweeps; s++)
    {
        gs_sweep_dia(&dia);
    }
    LIKWID_MARKER_STOP("DIA_sweep");
    auto t1DIA = std::chrono::high_resolution_clock::now();
    totalSweepTimeDIA += std::chrono::duration<double>(t1DIA - t0DIA).count();
    totalCallsDIA++;

    std::cout << "DIA_sweep: totalTime=" << totalSweepTimeDIA
              << "s calls=" << totalCallsDIA << std::endl;

    double res_dia = compute_residual(dia.psi, b, N);
    printf("DIA residual after %d sweeps: %.6e\n", num_sweeps, res_dia);
    free_dia(&dia);

    // RBDIA sweep loop with marker region and aggregate timing.
    RBDIAMatrix rbdia;
    assemble_rbdia(&rbdia, b, N);
    printf("\n--- RBDIA solver: %d sweeps ---\n", num_sweeps);

    static double totalSweepTimeRBDIA = 0.0;
    static int totalCallsRBDIA = 0;
    auto t0RBDIA = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("RBDIA_sweep");
    for (int s = 0; s < num_sweeps; s++)
    {
        gs_sweep_rbdia(&rbdia);
    }
    LIKWID_MARKER_STOP("RBDIA_sweep");
    auto t1RBDIA = std::chrono::high_resolution_clock::now();
    totalSweepTimeRBDIA += std::chrono::duration<double>(t1RBDIA - t0RBDIA).count();
    totalCallsRBDIA++;

    std::cout << "RBDIA_sweep: totalTime=" << totalSweepTimeRBDIA
              << "s calls=" << totalCallsRBDIA << std::endl;

    double res_rbdia = compute_residual(rbdia.psi, b, N);
    printf("RBDIA residual after %d sweeps: %.6e\n", num_sweeps, res_rbdia);
    free_rbdia(&rbdia);

    free(b);
    printf("\n=== Profiling complete ===\n");
}

/// @brief Run DIA sweeps until residual tolerance or sweep cap is reached.
/// @param N Grid resolution per dimension.
/// @param tol Target residual tolerance.
/// @param max_sweeps Maximum allowed sweep count.
void convergence_test_DIA(int N, double tol, int max_sweeps)
{
    int sweeps = 0;

    int nCells = N * N * N;

    // Allocate problem data and initialize DIA system.
    DIAMatrix diamat;
    double *bDia = (double *)malloc(sizeof(double) * nCells);
    double *uExactDia = (double *)malloc(sizeof(double) * nCells);

    compute_source(bDia, N);
    compute_exact(uExactDia, N);

    assemble_dia(&diamat, bDia, N);

    double residual = compute_residual(diamat.psi, bDia, N);
    printf("Initial residual = %.6e\n", residual);

    // Sweep until converged or capped; print residual periodically.
    while (residual > tol && sweeps < max_sweeps)
    {
        gs_sweep_dia(&diamat);
        sweeps++;

        if (sweeps % 50 == 0)
        {
            residual = compute_residual(diamat.psi, bDia, N);
            printf("Sweep %4d: residual = %.6e\n", sweeps, residual);
        }
    }

    // Report final accuracy against manufactured exact solution.
    double l2_error = compute_l2_error(diamat.psi, uExactDia, N);

    if (sweeps >= max_sweeps)
    {
        printf("Failed to converge in %d sweeps. Hit the cap.\n", sweeps);
        printf("N=%d did NOT converge, sweeps=%d, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }
    else
    {
        printf("N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }

    free_dia(&diamat);
    free(bDia);
    free(uExactDia);
}

/// @brief Run RBDIA sweeps until residual tolerance or sweep cap is reached.
/// @param N Grid resolution per dimension.
/// @param tol Target residual tolerance.
/// @param max_sweeps Maximum allowed sweep count.
void convergence_test_RBDIA(int N, double tol, int max_sweeps)
{
    int sweeps = 0;

    int nCells = N * N * N;

    // Allocate problem data and initialize RBDIA system.
    RBDIAMatrix rbdiamat;
    double *b = (double *)malloc(sizeof(double) * nCells);
    double *u_exact = (double *)malloc(sizeof(double) * nCells);

    compute_source(b, N);
    compute_exact(u_exact, N);

    assemble_rbdia(&rbdiamat, b, N);

    double residual = compute_residual(rbdiamat.psi, b, N);
    printf("Initial residual = %.6e\n", residual);

    // Sweep until converged or capped; print residual periodically.
    while (residual > tol && sweeps < max_sweeps)
    {
        gs_sweep_rbdia(&rbdiamat);
        sweeps++;

        if (sweeps % 50 == 0)
        {
            residual = compute_residual(rbdiamat.psi, b, N);
            printf("Sweep %4d: residual = %.6e\n", sweeps, residual);
        }
    }

    // Report final accuracy against manufactured exact solution.
    double l2_error = compute_l2_error(rbdiamat.psi, u_exact, N);

    if (sweeps >= max_sweeps)
    {
        printf("Failed to converge in %d sweeps. Hit the cap.\n", sweeps);
        printf("N=%d did NOT converge, sweeps=%d, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }
    else
    {
        printf("N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }

    free_rbdia(&rbdiamat);
    free(b);
    free(u_exact);
}

void convergence_test_LDU(int N, double tol, int max_sweeps)
{
    int sweeps = 0;

    int nCells = N * N * N;

    // Allocate problem data and initialize LDU system.
    LDUMatrix ldumat;
    double *b = (double *)malloc(sizeof(double) * nCells);
    double *u_exact = (double *)malloc(sizeof(double) * nCells);

    compute_source(b, N);
    compute_exact(u_exact, N);
    assemble_ldu(&ldumat, b, N);
    double residual = compute_residual(ldumat.psi, b, N);
    printf("Initial residual = %.6e\n", residual);

    // Sweep until converged or capped; print residual periodically.
    while (residual > tol && sweeps < max_sweeps)
    {
        gs_sweep_ldu(&ldumat);
        sweeps++;
        if (sweeps % 50 == 0)
        {
            residual = compute_residual(ldumat.psi, b, N);
            printf("Sweep %4d: residual = %.6e\n", sweeps, residual);
        }
    }

    // Always compute final residual so reporting is up-to-date.
    residual = compute_residual(ldumat.psi, b, N);

    // Report final accuracy against manufactured exact solution.
    double l2_error = compute_l2_error(ldumat.psi, u_exact, N);

    if (sweeps >= max_sweeps)
    {
        printf("Failed to converge in %d sweeps. Hit the cap.\n", sweeps);
        printf("N=%d did NOT converge, sweeps=%d, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }
    else
    {
        printf("N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error);
    }

    free_ldu(&ldumat);
    free(b);
    free(u_exact);
}

/// @brief Convenience wrapper to run both solver convergence tests sequentially.
/// @param N Grid resolution per dimension.
/// @param tol Target residual tolerance.
/// @param max_sweeps Maximum allowed sweep count.
void convergence_test_combined(int N, double tol, int max_sweeps)
{
    convergence_test_DIA(N, tol, max_sweeps);
    convergence_test_RBDIA(N, tol, max_sweeps);
}