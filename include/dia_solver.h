// No upper, lower, uAddr, ownerStart (this is how it differs from LDUMatrix)
// The off-diagonal coefficient is a single scalar (e.g. -1.0),
// not an array, because it's the same for every face
#include <stdio.h>
typedef struct {
    int N;              // grid dimension (square grid, so N×N cells)
    int nCells;         // number of cells (N*N)
    double *diag;       // diagonal coefficients, length N*N
    double *psi;        // solution vector, length N*N
    double *bPrime;     // working RHS, length N*N
    double *source;     // original RHS, length N*N
} DIAMatrix;

void assemble_dia(DIAMatrix *mat, double *source, int N);
void gs_sweep_dia(DIAMatrix *mat);
void gs_sweep_dia_trace(DIAMatrix *mat, FILE *fp);
void free_dia(DIAMatrix *mat);






// Function descriptions (brief)
// assemble_dia(DIAMatrix *mat, double *source, int N)
//     - allocate diag, psi, bPrime, source arrays (all length N*N)
//     - fill diag: 4.0 for every cell (Dirichlet BCs already baked in)
//     - copy source
//     - zero-initialise psi

// gs_sweep_dia(DIAMatrix *mat)
//     - reset bPrime = source
//     - interior sweep: i from 1 to N-2, j from 1 to N-2
//         pure offset arithmetic, no conditionals
//     - boundary sweeps: handle the four edges and four corners
//         (or do one combined boundary pass)

// free_dia(DIAMatrix *mat)
//     - free the 4 arrays
