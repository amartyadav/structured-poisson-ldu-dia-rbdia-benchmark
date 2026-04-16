#include <stdio.h>
typedef struct {
    int N;              // grid dimension (square grid, so N×N cells)
    int nCells;         // number of cells (N*N)
    double *diag;       // diagonal coefficients, length N*N
    double *psi;        // solution vector, length N*N
    double *bPrime;     // working RHS, length N*N
    double *source;     // original RHS, length N*N
} RBDIAMatrix;

void assemble_rbdia(RBDIAMatrix *mat, double *source, int N);
void gs_sweep_rbdia(RBDIAMatrix *mat);
void gs_sweep_rbdia_trace(RBDIAMatrix *mat, FILE *fp);
void free_rbdia(RBDIAMatrix *mat);
