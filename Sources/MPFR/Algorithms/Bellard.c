#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include "../../Headers/Sequential/Bellard_v1.h"


/*
 * Parallel Pi number calculation using the Bellard algorithm
 * Multiple threads can be used
 * The number of iterations is divided cyclically, 
 * so each thread calculates a part of Pi.  
 */
void Bellard_algorithm_OMP(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
    mpfr_t ONE; 

    mpfr_init_set_ui(ONE, 1, MPFR_RNDN); 

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, dep_a, dep_b, jump_dep_a, jump_dep_b, next_i;
        mpfr_t local_pi, dep_m, a, b, c, d, e, f, g, aux;

        thread_id = omp_get_thread_num();

        mpfr_init2(local_pi, precision_bits);               // private thread pi
        mpfr_set_ui(local_pi, 0, MPFR_RNDN);
        
        dep_a = thread_id * 4;
        dep_b = thread_id * 10;
        jump_dep_a = 4 * num_threads;
        jump_dep_b = 10 * num_threads;

        mpfr_init2(dep_m, precision_bits);
        mpfr_mul_2exp(dep_m, ONE, 10 * thread_id, MPFR_RNDN);
        mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);

        if(thread_id % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN);                   
        mpfr_inits2(precision_bits, a, b, c, d, e, f, g, aux, NULL);

        //First Phase -> Working on a local variable
        #pragma omp parallel for 
            for(i = thread_id; i < num_iterations; i+=num_threads){
                Bellard_iteration_v1(local_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
                // Update dependencies for next iteration:
                next_i = i + num_threads;
                mpfr_mul_2exp(dep_m, ONE, 10 * next_i, MPFR_RNDN);
                mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);
                if (next_i % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN); 
                dep_a += jump_dep_a;
                dep_b += jump_dep_b;  
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpfr_add(pi, pi, local_pi, MPFR_RNDN);

        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_pi, dep_m, a, b, c, d, e, f, g, aux, NULL);   
    }

    mpfr_div_ui(pi, pi, 64, MPFR_RNDN);
        
    //Clear memory
    mpfr_clear(ONE);
}

