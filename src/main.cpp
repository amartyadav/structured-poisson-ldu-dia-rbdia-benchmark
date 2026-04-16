#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "problem.hpp"
#include "ldu_solver.hpp"
#include "dia_solver.hpp"
#include "rb_dia_solver.hpp"
#include <chrono>


#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(a)
#define LIKWID_MARKER_STOP(a)
#define LIKWID_MARKER_CLOSE
#endif

void test_multiple_N(void);
void profile_sweeps(int N, int num_sweeps);
void convergence_test_DIA(int N, double tol, int max_sweeps);
void convergence_test_RBDIA(int N, double tol, int max_sweeps);
void convergence_test_combined(int N, double tol, int max_sweeps);


int main(int argc, char *argv[]) {
    if (argc > 1 && strcmp(argv[1], "profile") == 0)
    {
        int N = 100;
        int num_sweeps = 1000;
        if (argc > 2) N = atoi(argv[2]);
        if (argc > 3) num_sweeps = atoi(argv[3]);

        LIKWID_MARKER_INIT;
        LIKWID_MARKER_THREADINIT;

        profile_sweeps(N, num_sweeps);

        LIKWID_MARKER_CLOSE;
    }
    else if (argc > 1 && strcmp(argv[1], "convergence") == 0)
    {
        int N = 100;
        double tol = 1e-8;
        int max_sweeps = 10000;
        if (argc > 2) N = atoi(argv[2]);
        if (argc > 3) tol = atof(argv[3]);
        if (argc > 4) max_sweeps = atoi(argv[4]);

        convergence_test_combined(N, tol, max_sweeps);
    }
    else
    {
        test_multiple_N();
    }
    return 0;
}

void profile_sweeps(int N, int num_sweeps) {
    int nCells = N * N * N;
    int nFaces = 3 * N * N * (N - 1);

    printf("=== Profiling run: N=%d, sweeps=%d ===\n", N, num_sweeps);
    printf("Working set (DIA): ~%.1f MB\n",
           (4.0 * nCells * sizeof(double)) / (1024.0 * 1024.0));
    // printf("Working set (LDU): ~%.1f MB\n",
    //        (4.0 * nCells * sizeof(double)
    //         + 3.0 * nFaces * sizeof(double)
    //         + nFaces * sizeof(int)
    //         + (nCells + 1) * sizeof(int))
    //        / (1024.0 * 1024.0));

    double *b = (double*)malloc(sizeof(double) * nCells);
    compute_source(b, N);

    // --- LDU solver ---
    // LDUMatrix ldu;
    // assemble_ldu(&ldu, b, N);
    // printf("\n--- LDU solver: %d sweeps ---\n", num_sweeps);

    // LIKWID_MARKER_START("LDU_sweep");
    // for (int s = 0; s < num_sweeps; s++)
    // {
    //     static double totalSweepTimeLDU = 0.0;
    //     static int totalCallsLDU = 0;
    //     auto t0LDU = std::chrono::high_resolution_clock::now();
    //     gs_sweep_ldu(&ldu);
    //     auto t1LDU = std::chrono::high_resolution_clock::now();
    //     totalSweepTimeLDU += std::chrono::duration<double>(t1LDU - t0LDU).count();
    //     totalCallsLDU++;
    //     if (totalCallsLDU % 1000 == 0)
    //     {
    //         std::cout << "LDU_sweep: totalTime=" << totalSweepTimeLDU
    //             << "s calls=" << totalCallsLDU << std::endl;
    //     }
    // }
    // LIKWID_MARKER_STOP("LDU_sweep");


    // double res_ldu = compute_residual(ldu.psi, b, N);
    // printf("LDU residual after %d sweeps: %.6e\n", num_sweeps, res_ldu);
    // free_ldu(&ldu);

    // --- DIA solver ---
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


    // --- RBDIA solver ---
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

void test_multiple_N(void) {
    int N_s[] = {4, 6, 8, 12, 16, 20};
    int num_sizes = 6;

    printf("=== LDU Solver Runs ===\n\n");
    for (int i = 0; i < num_sizes; i++)
    {
        int N = N_s[i];
        int nCells = N * N * N;
        printf("=== Starting LDU GS for N=%d ===\n", N);
        LDUMatrix mat;
        double *b = (double *)malloc(sizeof(double) * nCells);
        double *u_exact = (double *)malloc(sizeof(double) * nCells);
        compute_source(b, N);
        compute_exact(u_exact, N);
        assemble_ldu(&mat, b, N);
        printf("diag[0]=%.1f  diag[nCells/2]=%.1f  diag[nCells-1]=%.1f\n",
               mat.diag[0], mat.diag[nCells / 2], mat.diag[nCells - 1]);
        printf("nFaces expected=%d  ownerStart[nCells]=%d\n",
               mat.nFaces, mat.ownerStart[mat.nCells]);
        for (int sweep = 0; sweep < 1000; sweep++)
        {
            static double totalSweepTimeLDU = 0.0;
            static int totalCallsLDU = 0;
            auto t0LDU = std::chrono::high_resolution_clock::now();
            gs_sweep_ldu(&mat);
            auto t1LDU = std::chrono::high_resolution_clock::now();
            totalSweepTimeLDU += std::chrono::duration<double>(t1LDU - t0LDU).count();
            totalCallsLDU++;
            if (totalCallsLDU % 1000 == 0)
            {
                std::cout << "LDU_sweep: totalTime=" << totalSweepTimeLDU
                    << "s calls=" << totalCallsLDU << std::endl;
            }
            if (sweep % 500 == 0)
            {
                double res = compute_residual(mat.psi, b, N);
                printf("Sweep %4d: residual = %.6e\n", sweep, res);
                if (res < 1e-14)
                {
                    break;
                }
            }
        }
        double l2 = compute_l2_error(mat.psi, u_exact, N);
        printf("L2 error = %.6e\n", l2);
        free_ldu(&mat);
        free(b);
        free(u_exact);
        printf("=== === ===\n");
    }

    printf("\n=== DIA Solver Runs ===\n\n");
    for (int i = 0; i < num_sizes; i++)
    {
        int N = N_s[i];
        int nCells = N * N * N;
        printf("=== Starting DIA GS for N=%d ===\n", N);
        DIAMatrix dmat;
        double *b = (double *)malloc(sizeof(double) * nCells);
        double *u_exact = (double *)malloc(sizeof(double) * nCells);
        compute_source(b, N);
        compute_exact(u_exact, N);
        assemble_dia(&dmat, b, N);
        printf("diag[0]=%.1f  diag[nCells/2]=%.1f  diag[nCells-1]=%.1f\n",
               dmat.diag[0], dmat.diag[nCells / 2], dmat.diag[nCells - 1]);
        for (int sweep = 0; sweep < 10000; sweep++)
        {
            static double totalSweepTimeDIA = 0.0;
            static int totalCallsDIA = 0;
            auto t0DIA = std::chrono::high_resolution_clock::now();
            gs_sweep_dia(&dmat);
            auto t1DIA = std::chrono::high_resolution_clock::now();
            totalSweepTimeDIA += std::chrono::duration<double>(t1DIA - t0DIA).count();
            totalCallsDIA++;
            if (totalCallsDIA % 1000 == 0)
            {
                std::cout << "DIA_sweep: totalTime=" << totalSweepTimeDIA
                    << "s calls=" << totalCallsDIA << std::endl;
            }
            if (sweep % 500 == 0)
            {
                double res = compute_residual(dmat.psi, b, N);
                printf("Sweep %4d: residual = %.6e\n", sweep, res);
                if (res < 1e-14)
                {
                    break;
                }
            }
        }
        double l2 = compute_l2_error(dmat.psi, u_exact, N);
        printf("L2 error = %.6e\n", l2);
        free_dia(&dmat);
        free(b);
        free(u_exact);
        printf("=== === ===\n");
    }

    printf("\n=== RBDIA Solver Runs ===\n\n");
    for (int i = 0; i < num_sizes; i++)
    {
        int N = N_s[i];
        int nCells = N * N * N;
        printf("=== Starting RBDIA GS for N=%d ===\n", N);
        RBDIAMatrix dmat;
        double *b = (double *)malloc(sizeof(double) * nCells);
        double *u_exact = (double *)malloc(sizeof(double) * nCells);
        compute_source(b, N);
        compute_exact(u_exact, N);
        assemble_rbdia(&dmat, b, N);
        printf("diag[0]=%.1f  diag[nCells/2]=%.1f  diag[nCells-1]=%.1f\n",
               dmat.diag[0], dmat.diag[nCells / 2], dmat.diag[nCells - 1]);
        for (int sweep = 0; sweep < 10000; sweep++)
        {
            static double totalSweepTimeRBDIA = 0.0;
            static int totalCallsRBDIA = 0;
            auto t0RBDIA = std::chrono::high_resolution_clock::now();
            gs_sweep_rbdia(&dmat);
            auto t1RBDIA = std::chrono::high_resolution_clock::now();
            totalSweepTimeRBDIA += std::chrono::duration<double>(t1RBDIA - t0RBDIA).count();
            totalCallsRBDIA++;
            if (totalCallsRBDIA % 1000 == 0)
            {
                std::cout << "RBDIA_sweep: totalTime=" << totalSweepTimeRBDIA
                    << "s calls=" << totalCallsRBDIA << std::endl;
            }
            if (sweep % 500 == 0)
            {
                double res = compute_residual(dmat.psi, b, N);
                printf("Sweep %4d: residual = %.6e\n", sweep, res);
                if (res < 1e-14)
                {
                    break;
                }
            }
        }
        double l2 = compute_l2_error(dmat.psi, u_exact, N);
        printf("L2 error = %.6e\n", l2);
        free_rbdia(&dmat);
        free(b);
        free(u_exact);
        printf("=== === ===\n");
    }
}

void convergence_test_DIA(int N, double tol, int max_sweeps)
{
    int sweeps = 0;

    int nCells = N * N * N;

    // DIA Test
    DIAMatrix diamat;
    double *bDia = (double *)malloc(sizeof(double) * nCells);
    double *uExactDia = (double *)malloc(sizeof(double) * nCells);


    compute_source(bDia, N);
    compute_exact(uExactDia, N);

    assemble_dia(&diamat, bDia, N);

    double residual = compute_residual(diamat.psi, bDia, N);
    printf("Initial residual = %.6e\n", residual);

    while(residual > tol && sweeps < max_sweeps)
    {
        gs_sweep_dia(&diamat);
        sweeps++;

        if(sweeps % 50 == 0)
        {
            residual = compute_residual(diamat.psi, bDia, N);
            printf("Sweep %4d: residual = %.6e\n", sweeps, residual);
        }
    }

    if (sweeps >= max_sweeps)
    {
        printf("Failed to converge in %d sweeps. Hit the cap.\n", sweeps);
    }

    // computing L2 error
    double l2_error = compute_l2_error(diamat.psi, uExactDia, N);

    // "N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e"
    if(residual < tol) { printf("N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error); }


    free_dia(&diamat);
    free(bDia);
    free(uExactDia);
}

void convergence_test_RBDIA(int N, double tol, int max_sweeps)
{
    int sweeps = 0;

    int nCells = N * N * N;

    // RBDIA Test
    RBDIAMatrix rbdiamat;
    double *b = (double *)malloc(sizeof(double) * nCells);
    double *u_exact = (double *)malloc(sizeof(double) * nCells);


    compute_source(b, N);
    compute_exact(u_exact, N);

    assemble_rbdia(&rbdiamat, b, N);

    double residual = compute_residual(rbdiamat.psi, b, N);
    printf("Initial residual = %.6e\n", residual);

    while(residual > tol && sweeps < max_sweeps)
    {
        gs_sweep_rbdia(&rbdiamat);
        sweeps++;

        if(sweeps % 50 == 0)
        {
            residual = compute_residual(rbdiamat.psi, b, N);
            printf("Sweep %4d: residual = %.6e\n", sweeps, residual);
        }
    }

    if (sweeps >= max_sweeps)
    {
        printf("Failed to converge in %d sweeps. Hit the cap.\n", sweeps);
    }

    // computing L2 error
    double l2_error = compute_l2_error(rbdiamat.psi, u_exact, N);

    // "N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e"
    if(residual < tol) { printf("N=%d converged in %d sweeps, final residual %.6e, L2 error %.6e\n", N, sweeps, residual, l2_error); }

    free_rbdia(&rbdiamat);
    free(b);
    free(u_exact);
}

void convergence_test_combined(int N, double tol, int max_sweeps)
{
    convergence_test_DIA(N, tol, max_sweeps);
    convergence_test_RBDIA(N, tol, max_sweeps);
}
