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
 * This version uses mpz to make faster the computations related with integers      *
 * This version uses a blocks distribution                                          *
 * It allows to compute pi using multiple threads                                   *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula:                                                              *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula dependencies:                                                 *
 ************************************************************************************/


/*
 * An iteration of Chudnovsky formula using the binary splitting algorithm
 */
void gmp_chudnovsky_craig_wood_expression_iteration(int n, mpf_t sum_a, mpf_t sum_b, mpf_t a_n, mpf_t b_n, mpf_t a_n_divisor){
    // Computing a_n
    mpf_mul_ui(a_n, a_n, - (((6 * n) - 1) *((2 * n) - 1) * ((6 * n) - 1)));
    mpf_set_ui(a_n_divisor, n);
    mpf_pow_ui(a_n_divisor, a_n_divisor, 3);
    mpf_mul_ui(a_n_divisor, a_n_divisor, C);
    mpf_div(a_n, a_n, a_n_divisor);

    // Computing b_n
    mpf_mul_ui(b_n, a_n, n);

    // Computing sum_a and sum_b
    mpf_add(sum_a, sum_a, a_n);
    mpf_add(sum_b, sum_b, b_n);

    gmp_printf(" sum_a ----- %.Ff \n", sum_a);
    gmp_printf(" sum_b ----- %.Ff \n", sum_b);
}

void gmp_chudnovsky_craig_wood_expression_blocks_algorithm(mpf_t pi, int num_iterations, int num_threads){
    mpf_t e;
    mpf_init_set_ui(e, E);

    int i;
    mpf_t sum_a, sum_b, a_n, b_n, a_n_divisor;

    mpf_inits(sum_a, sum_b, a_n, b_n, a_n_divisor, NULL);

    mpf_set_ui(a_n, 1);
    mpf_set_ui(sum_a, 1);
    mpf_set_ui(sum_b, 0);

    for(i = 1; i < num_iterations; i++){
        gmp_chudnovsky_craig_wood_expression_iteration(i, sum_a, sum_b, a_n, b_n, a_n_divisor);
    }
    

    mpf_mul_ui(sum_a, sum_a, A);
    mpf_mul_ui(sum_b, sum_b, B);
    mpf_add(pi, sum_a, sum_b);
    
    mpf_sqrt(e, e);
    mpf_mul_ui(e, e, D);
    mpf_div(pi, e, pi);    
        
    //Clear thread memory
    mpf_clears(sum_a, sum_b, a_n, b_n, a_n_divisor, NULL);  
   
    //Clear memory
    mpf_clear(e);
}