typedef struct
{
    double *diag;
    double *upper;
    double *lower;
    int *uAddr;
    int *ownerStart;
    double *psi;
    double *bPrime;
    double *source;
    int nCells;
    int nFaces;
} LDUMatrix;

void assemble_ldu(LDUMatrix *mat, double *source, int N);
void gs_sweep_ldu(LDUMatrix *mat);
void free_ldu(LDUMatrix *mat);
