#ifndef CHUDNOVSKY_GMP
#define CHUDNOVSKY_GMP

void chudnovsky_algorithm_gmp(mpf_t pi, int num_iterations, int num_threads);
void init_dep_a(mpf_t dep_a, int block_start);

#endif

