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
void Bellard_algorithm_v1_OMP(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
    mpfr_t jump; 

    mpfr_init2(jump, precision_bits);
    mpfr_set_ui(jump, 1, MPFR_RNDN); 
    mpfr_div_ui(jump, jump, 1024, MPFR_RNDN);
    mpfr_pow_ui(jump, jump, num_threads, MPFR_RNDN);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, dep_a, dep_b, jump_dep_a, jump_dep_b;
        mpfr_t local_pi, dep_m, a, b, c, d, e, f, g, aux;

        thread_id = omp_get_thread_num();

        mpfr_init2(local_pi, precision_bits);               // private thread pi
        mpfr_set_ui(local_pi, 0, MPFR_RNDN);
        dep_a = thread_id * 4;
        dep_b = thread_id * 10;
        jump_dep_a = 4 * num_threads;
        jump_dep_b = 10 * num_threads;
        mpfr_init2(dep_m, precision_bits);
        mpfr_set_ui(dep_m, 1, MPFR_RNDN);
        mpfr_div_ui(dep_m, dep_m, 1024, MPFR_RNDN);
        mpfr_pow_ui(dep_m, dep_m, thread_id, MPFR_RNDN);        // dep_m = ((-1)^n)/1024)
        if(thread_id % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN);                   
        mpfr_inits2(precision_bits, a, b, c, d, e, f, g, aux, NULL);

        //First Phase -> Working on a local variable
        if(num_threads % 2 != 0){
            #pragma omp parallel for 
                for(i = thread_id; i < num_iterations; i+=num_threads){
                    Bellard_iteration_v1(local_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
                    // Update dependencies for next iteration:
                    mpfr_mul(dep_m, dep_m, jump, MPFR_RNDN); 
                    mpfr_neg(dep_m, dep_m, MPFR_RNDN); 
                    dep_a += jump_dep_a;
                    dep_b += jump_dep_b;  
                }
        } else {
            #pragma omp parallel for
                for(i = thread_id; i < num_iterations; i+=num_threads){
                    Bellard_iteration_v1(local_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
                    // Update dependencies for next iteration:
                    mpfr_mul(dep_m, dep_m, jump, MPFR_RNDN);    
                    dep_a += jump_dep_a;
                    dep_b += jump_dep_b;  
                }
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
    mpfr_clear(jump);
}

