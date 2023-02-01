#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include "bellard_recursive_power_cyclic.h"

double gettimeofday();

/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Bellard formula implementation                                                   *
 * It implements a single-threaded method and another that can use multiple threads *
 * It uses a cyclic distribution                                                    *
 *                                                                                  *
 ************************************************************************************
 * Bellard formula:                                                                 *
 *                 (-1)^n     32     1      256     64       4       4       1      *
 * 2^6 * pi = SUM( ------ [- ---- - ---- + ----- - ----- - ----- - ----- + -----])  *
 *                 1024^n    4n+1   4n+3   10n+1   10n+3   10n+5   10n+7   10n+9    *
 *                                                                                  *
 * Formula quotients are coded as:                                                  *
 *             32          1           256          64                              *
 *        a = ----,   b = ----,   c = -----,   d = -----,                           *
 *            4n+1        4n+3        10n+1        10n+3                            *
 *                                                                                  *
 *              4            4            1         (-1)^n                          *
 *        e = -----,   f = -----,   g = -----,   m = -----,                         *
 *            10n+5        10n+7        10n+9        2^10n                          *
 *                                                                                  *
 ************************************************************************************
 * Bellard formula dependencies:                                                    *
 *                           1            1                                         *
 *              dep_m(n) = ------ = -----------------                               *
 *                         1024^n   1024^(n-1) * 1024                               *
 *                                                                                  *
 *              dep_a(n) = 4n  = dep_a(n-1) + 4                                     *
 *                                                                                  *
 *              dep_b(n) = 10n = dep_a(n-1) + 10                                    *
 *                                                                                  *
 ************************************************************************************/


void mpfr_bellard_bit_shift_power_cyclic_algorithm(mpfr_t pi, int num_iterations, int num_threads){
    mpfr_t ONE; 
    mpfr_init_set_ui(ONE, 1, MPFR_RNDN); 

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, dep_a, dep_b, jump_dep_a, jump_dep_b, next_i;
        mpfr_t local_pi, dep_m, a, b, c, d, e, f, g, aux;

        thread_id = omp_get_thread_num();
        mpfr_init_set_ui(local_pi, 0, MPFR_RNDN);
        dep_a = thread_id * 4;
        dep_b = thread_id * 10;
        jump_dep_a = 4 * num_threads;
        jump_dep_b = 10 * num_threads;
        mpfr_init(dep_m);
        mpfr_mul_2exp(dep_m, ONE, 10 * thread_id, MPFR_RNDN);
        mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);
        if(thread_id % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN);                   
        mpfr_inits(a, b, c, d, e, f, g, aux, NULL);

        double execution_time;
        struct timeval t1, t2;
        //First Phase -> Working on a local variable
        for(i = thread_id; i < num_iterations; i+=num_threads){
            gettimeofday(&t1, NULL);
            mpfr_bellard_iteration(local_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
            // Update dependencies for next iteration:
            next_i = i + num_threads;
            mpfr_mul_2exp(dep_m, ONE, 10 * next_i, MPFR_RNDN);
            mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);
            if (next_i % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN); 
            dep_a += jump_dep_a;
            dep_b += jump_dep_b;  
            gettimeofday(&t2, NULL);
            execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e3; 
            printf("%f\n", execution_time);
        }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpfr_add(pi, pi, local_pi, MPFR_RNDN);

        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_pi, dep_m, a, b, c, d, e, f, g, aux, NULL);   
    }

    mpfr_div_2exp(pi, pi, 6, MPFR_RNDN); // pi = pi / 2⁶
        
    //Clear memory
    mpfr_clear(ONE);
}

