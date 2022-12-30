#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
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
 * An iteration of Chudnovsky formula
 */
void mpfr_chudnovsky_craig_wood_expression_iteration(int n, mpfr_t sum_a, mpfr_t sum_b, mpfr_t a_n, mpfr_t b_n, mpfr_t a_n_divisor, 
                                                    mpfr_t factor_a, mpfr_t factor_b, mpfr_t factor_c){
    // Computing a_n
    // factor_a = (6n -5)
    mpfr_set_ui(factor_a, n, MPFR_RNDN);
    mpfr_mul_ui(factor_a, factor_a, 6, MPFR_RNDN);
    mpfr_sub_ui(factor_a, factor_a, 5, MPFR_RNDN);
    
    // factor_a = (2n -1)
    mpfr_set_ui(factor_b, n, MPFR_RNDN);
    mpfr_mul_ui(factor_b, factor_b, 2, MPFR_RNDN);
    mpfr_sub_ui(factor_b, factor_b, 1, MPFR_RNDN);

    // factor_a = (6n -1)
    mpfr_set_ui(factor_c, n, MPFR_RNDN);
    mpfr_mul_ui(factor_c, factor_c, 6, MPFR_RNDN);
    mpfr_sub_ui(factor_c, factor_c, 1, MPFR_RNDN);

    // factor_a = -((6n -5) * (2n - 1) * (6n - 1))
    mpfr_mul(factor_a, factor_a, factor_b, MPFR_RNDN);
    mpfr_mul(factor_a, factor_a, factor_c, MPFR_RNDN);
    mpfr_neg(factor_a, factor_a, MPFR_RNDN);

    // a_n = (factor_a * a_n_prev) / (n^3 * C)
    mpfr_mul(a_n, a_n, factor_a, MPFR_RNDN);
    mpfr_set_ui(a_n_divisor, n, MPFR_RNDN);
    mpfr_pow_ui(a_n_divisor, a_n_divisor, 3, MPFR_RNDN);
    mpfr_mul_ui(a_n_divisor, a_n_divisor, C, MPFR_RNDN);
    mpfr_div(a_n, a_n, a_n_divisor, MPFR_RNDN);

    // Computing b_n
    // b_n = a_n * n
    mpfr_mul_ui(b_n, a_n, n, MPFR_RNDN);

    // Computing sum_a and sum_b
    mpfr_add(sum_a, sum_a, a_n, MPFR_RNDN);
    mpfr_add(sum_b, sum_b, b_n, MPFR_RNDN);

}


void mpfr_chudnovsky_craig_wood_expression_blocks_algorithm(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
    mpfr_t sum_a, sum_b, e;

    mpfr_inits2(precision_bits, sum_a, sum_b, e, NULL);
    mpfr_set_ui(e, E, MPFR_RNDN);
    mpfr_set_ui(sum_a, 1, MPFR_RNDN);
    mpfr_set_ui(sum_b, 0, MPFR_RNDN);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {   
        int thread_id, i, block_start, block_end;
        mpfr_t a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c;

        thread_id = omp_get_thread_num();
        mpfr_inits2(precision_bits, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c, NULL);
        mpfr_set_ui(a_n, 1, MPFR_RNDN);
        

        for(i = 1; i < num_iterations; i++){
            mpfr_chudnovsky_craig_wood_expression_iteration(i, sum_a, sum_b, a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c);
        }

        //Clear thread memory
        mpfr_clears(a_n, b_n, a_n_divisor, factor_a, factor_b, factor_c, NULL);   
    }

    mpfr_mul_ui(sum_a, sum_a, A, MPFR_RNDN);
    mpfr_mul_ui(sum_b, sum_b, B, MPFR_RNDN);
    mpfr_add(pi, sum_a, sum_b, MPFR_RNDN);
    
    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
   
    //Clear memory
    mpfr_clears(sum_a, sum_b, e, NULL);
}