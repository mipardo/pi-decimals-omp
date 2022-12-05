#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "chudnovsky_blocks_with_all_factorials.h"

#define A 13591409
#define B 545140134
#define C 640320
#define D 426880
#define E 10005


/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Chudnovsky formula implementation                                                *
 * This version does not computes all the factorials (the expression is simplified) *
 * This version uses a block distribution                                           *
 * It allows to compute pi using multiple threads                                   *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula:                                                              *
 *     426880 sqrt(10005)                 (6n)! (545140134n + 13591409)             *
 *    --------------------  = SUMMATORY( ----------------------------- ),  n >=0    *
 *            pi                            (n!)^3 (3n)! (-640320)^3n               *
 *                                                                                  *
 * Some operands of the formula are coded as:                                       *
 *      dep_a_dividend = (6n)!                                                      *
 *      dep_a_divisor  = (n!)^3 (3n)!                                               *
 *      e              = 426880 sqrt(10005)                                         *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula dependencies:                                                 *
 *                     (6n)!         (12n + 10)(12n + 6)(12n + 2)                   *
 *      dep_a(n) = --------------- = ---------------------------- * dep_a(n-1)      *
 *                 ((n!)^3 (3n)!)              (n + 1)^3                            *
 *                                                                                  *
 *      dep_b(n) = (-640320)^3n = (-640320)^3(n-1) * (-640320)^3)                   *
 *                                                                                  *
 *      dep_c(n) = (545140134n + 13591409) = dep_c(n - 1) + 545140134               *
 *                                                                                  *
 ************************************************************************************/

/*
 * An iteration of Chudnovsky formula
 */
void chudnovsky_iteration_gmp(mpf_t pi, int n, mpf_t dep_a, mpf_t dep_b, mpf_t dep_c, mpf_t aux){
    mpf_mul(aux, dep_a, dep_c);
    mpf_div(aux, aux, dep_b);
    
    mpf_add(pi, pi, aux);
}

/*
 * This method is used by ParallelChudnovskyAlgorithm threads
 * for computing the first value of dep_a
 */
void init_dep_a_gmp(mpf_t dep_a, int block_start){
    mpz_t factorial_n, dividend, divisor;
    mpf_t float_dividend, float_divisor;
    mpz_inits(factorial_n, dividend, divisor, NULL);
    mpf_inits(float_dividend, float_divisor, NULL);

    mpz_fac_ui(factorial_n, block_start);
    mpz_fac_ui(divisor, 3 * block_start);
    mpz_fac_ui(dividend, 6 * block_start);

    mpz_pow_ui(factorial_n, factorial_n, 3);
    mpz_mul(divisor, divisor, factorial_n);

    mpf_set_z(float_dividend, dividend);
    mpf_set_z(float_divisor, divisor);

    mpf_div(dep_a, float_dividend, float_divisor);

    mpz_clears(factorial_n, dividend, divisor, NULL);
    mpf_clears(float_dividend, float_divisor, NULL);
}

void chudnovsky_blocks_with_simplified_expression_algorithm_gmp(mpf_t pi, int num_iterations, int num_threads){
    mpf_t e, c;
    int block_size;
    
    block_size = (num_iterations + num_threads - 1) / num_threads;
    mpf_init_set_ui(e, E);
    mpf_init_set_ui(c, C);
    mpf_neg(c, c);
    mpf_pow_ui(c, c, 3);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {   
        int thread_id, i, block_start, block_end, factor_a;
        mpf_t local_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux;

        thread_id = omp_get_thread_num();
        block_start = thread_id * block_size;
        block_end = block_start + block_size;
        if (block_end > num_iterations) block_end = num_iterations;
        
        mpf_init_set_ui(local_pi, 0);    // private thread pi
        mpf_inits(dep_a, dep_b, dep_a_dividend, dep_a_divisor, aux, NULL);
        init_dep_a_gmp(dep_a, block_start);
        mpf_pow_ui(dep_b, c, block_start);
        mpf_init_set_ui(dep_c, B);
        mpf_mul_ui(dep_c, dep_c, block_start);
        mpf_add_ui(dep_c, dep_c, A);
        factor_a = 12 * block_start;

        //First Phase -> Working on a local variable        
        for(i = block_start; i < block_end; i++){
            chudnovsky_iteration_gmp(local_pi, i, dep_a, dep_b, dep_c, aux);
            //Update dep_a:
            mpf_set_ui(dep_a_dividend, factor_a + 10);
            mpf_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 6);
            mpf_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 2);
            mpf_mul(dep_a_dividend, dep_a_dividend, dep_a);

            mpf_set_ui(dep_a_divisor, i + 1);
            mpf_pow_ui(dep_a_divisor, dep_a_divisor ,3);
            mpf_div(dep_a, dep_a_dividend, dep_a_divisor);
            factor_a += 12;

            //Update dep_b:
            mpf_mul(dep_b, dep_b, c);

            //Update dep_c:
            mpf_add_ui(dep_c, dep_c, B);
        }

        //Second Phase -> Accumulate the result in the global variable 
        #pragma omp critical
        mpf_add(pi, pi, local_pi);
        
        //Clear thread memory
        mpf_clears(local_pi, dep_a, dep_b, dep_c, dep_a_dividend, dep_a_divisor, aux, NULL);   
    }

    mpf_sqrt(e, e);
    mpf_mul_ui(e, e, D);
    mpf_div(pi, e, pi);    
    
    //Clear memory
    mpf_clears(c, e, NULL);
}