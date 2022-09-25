#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include "../../Headers/Sequential/BBP.h"

#define QUOTIENT 0.0625


/*
 * Parallel Pi number calculation using the BBP algorithm
 * Multiple threads can be used
 * The number of iterations is divided in blocks, 
 * so each thread calculates a part of Pi.  
 */

void BBP_algorithm_OMP(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
    mpfr_t quotient; 

    mpfr_init_set_d(quotient, QUOTIENT, MPFR_RNDN);         // quotient = (1 / 16)   

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, block_size, block_start, block_end;
        mpfr_t local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux;

        thread_id = omp_get_thread_num();
        block_size = (num_iterations + num_threads - 1) / num_threads;
        block_start = thread_id * block_size;
        block_end = block_start + block_size;
        if (block_end > num_iterations) block_end = num_iterations;
        
        mpfr_init2(local_pi, precision_bits);               // private thread pi
        mpfr_set_ui(local_pi, 0, MPFR_RNDN);
        mpfr_init2(dep_m, precision_bits);
        mpfr_pow_ui(dep_m, quotient, block_start, MPFR_RNDN);    // m = (1/16)^n                  
        mpfr_inits2(precision_bits, quot_a, quot_b, quot_c, quot_d, aux, NULL);
        

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = block_start; i < block_end; i++){
                BBP_iteration(local_pi, i, dep_m, quot_a, quot_b, quot_c, quot_d, aux);
                // Update dependencies:  
                mpfr_mul(dep_m, dep_m, quotient, MPFR_RNDN);
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpfr_add(pi, pi, local_pi, MPFR_RNDN);


        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux, NULL);   
    }
        
    //Clear memory
    mpfr_clear(quotient);
}