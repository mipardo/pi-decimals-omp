#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>

#define A 13591409
#define B 545140134
#define C 10939058860032000
#define D 426880
#define E 10005


/************************************************************************************
 * Miguel Pardo Navarro. 25/12/2022                                                 *
 * Chudnovsky formula implementation                                                *
 * This version use the Craig Wood simplified mathematical expression               * 
 * See https://www.craig-wood.com/nick/articles/pi-chudnovsky                       * 
 * The iterations of the series can not be paralelized                              *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula:                                                              *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula dependencies:                                                 *
 ************************************************************************************/

double gettimeofday();

/*
 * An iteration of Chudnovsky formula using the binary splitting algorithm
 */
void gmp_chudnovsky_craig_wood_expression_iteration(int n, mpf_t sum_a, mpf_t sum_b, mpf_t a_n, mpf_t b_n, mpf_t a_n_divisor, 
                                                    mpf_t factor_a, mpf_t factor_b, mpf_t factor_c){
    // Computing a_n
    // factor_a = (6n -5)
    mpf_set_ui(factor_a, n);
    mpf_mul_ui(factor_a, factor_a, 6);
    mpf_sub_ui(factor_a, factor_a, 5);
    
    // factor_a = (2n -1)
    mpf_set_ui(factor_b, n);
    mpf_mul_ui(factor_b, factor_b, 2);
    mpf_sub_ui(factor_b, factor_b, 1);

    // factor_a = (6n -1)
    mpf_set_ui(factor_c, n);
    mpf_mul_ui(factor_c, factor_c, 6);
    mpf_sub_ui(factor_c, factor_c, 1);

    // factor_a = -((6n -5) * (2n - 1) * (6n - 1))
    mpf_mul(factor_a, factor_a, factor_b);
    mpf_mul(factor_a, factor_a, factor_c);
    mpf_neg(factor_a, factor_a);

    // a_n = (factor_a * a_n_prev) / (n^3 * C)
    mpf_mul(a_n, a_n, factor_a);
    mpf_set_ui(a_n_divisor, n);
    mpf_pow_ui(a_n_divisor, a_n_divisor, 3);
    mpf_mul_ui(a_n_divisor, a_n_divisor, C);
    mpf_div(a_n, a_n, a_n_divisor);

    // Computing b_n
    // b_n = a_n * n
    mpf_mul_ui(b_n, a_n, n);

    // Computing sum_a and sum_b
    mpf_add(sum_a, sum_a, a_n);
    mpf_add(sum_b, sum_b, b_n);
}


void gmp_chudnovsky_craig_wood_expression_algorithm(mpf_t pi, int num_iterations, int num_threads){
    int i;
    mpf_t sum_a, sum_b, e, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c;

    mpf_inits(sum_a, sum_b, e, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c, NULL);
    mpf_set_ui(a_n, 1);
    mpf_set_ui(e, E);
    mpf_set_ui(sum_a, 1);
    mpf_set_ui(sum_b, 0);

    double execution_time;
    struct timeval t1, t2;
    for(i = 1; i < num_iterations; i++){
        gettimeofday(&t1, NULL);
        gmp_chudnovsky_craig_wood_expression_iteration(i, sum_a, sum_b, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c);
        gettimeofday(&t2, NULL);
        execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e3; 
        printf("%f\n", execution_time);
    }

    mpf_mul_ui(sum_a, sum_a, A);
    mpf_mul_ui(sum_b, sum_b, B);
    mpf_add(pi, sum_a, sum_b);
    
    mpf_sqrt(e, e);
    mpf_mul_ui(e, e, D);
    mpf_div(pi, e, pi);    
   
    //Clear memory
    mpf_clears(sum_a, sum_b, e, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c, NULL);
}