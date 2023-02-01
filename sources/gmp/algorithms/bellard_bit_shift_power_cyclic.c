#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "bellard_recursive_power_cyclic.h"


/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Bellard formula implementation                                                   *
 * It allows to compute pi using multiple threads                                   *
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

double gettimeofday();

void gmp_bellard_bit_shift_power_cyclic_algorithm(mpf_t pi, int num_iterations, int num_threads){
    mpf_t ONE;
    mpf_init_set_ui(ONE, 1);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, dep_a, dep_b, jump_dep_a, jump_dep_b, next_i;
        mpf_t local_pi, dep_m, a, b, c, d, e, f, g, aux;

        thread_id = omp_get_thread_num();
        mpf_init_set_ui(local_pi, 0);       // private thread pi
        dep_a = thread_id * 4;
        dep_b = thread_id * 10;
        jump_dep_a = 4 * num_threads;
        jump_dep_b = 10 * num_threads;
        mpf_init(dep_m);
        mpf_mul_2exp(dep_m, ONE, 10 * thread_id);
        mpf_div(dep_m, ONE, dep_m);
        if(thread_id % 2 != 0) mpf_neg(dep_m, dep_m);                   
        mpf_inits(a, b, c, d, e, f, g, aux, NULL);

        double execution_time;
        struct timeval t1, t2;
        //First Phase -> Working on a local variable
        for(i = thread_id; i < num_iterations; i+=num_threads){
            gettimeofday(&t1, NULL);
            gmp_bellard_iteration(local_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
            // Update dependencies for next iteration:
            next_i = i + num_threads;
            mpf_mul_2exp(dep_m, ONE, 10 * next_i);
            mpf_div(dep_m, ONE, dep_m);
            if (next_i % 2 != 0) mpf_neg(dep_m, dep_m); 
            dep_a += jump_dep_a;
            dep_b += jump_dep_b;
            gettimeofday(&t2, NULL);
            execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e3; 
            printf("%f\n", execution_time);
        }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpf_add(pi, pi, local_pi);

        //Clear thread memory
        mpf_clears(local_pi, dep_m, a, b, c, d, e, f, g, aux, NULL);   
    }

    mpf_div_2exp(pi, pi, 6); // pi = pi / 2⁶

    mpf_clear(ONE);        
}

