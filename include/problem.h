void compute_source(double *b, int N);

void compute_exact(double *u_exact, int N);

double compute_l2_error(double *u, double *u_exact, int N);

double compute_residual(double *u, double *b, int N);
