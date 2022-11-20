#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "bbp_cyclic.h"


#define QUOTIENT 0.0625

/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Last version of Bailey Borwein Plouffe formula implementation                    *
 * It allows to compute pi using multiple threads                                   *
 * It uses a block distribution                                                     *
 *                                                                                  *
 ************************************************************************************
 * Bailey Borwein Plouffe formula:                                                  *
 *                      1        4          2        1       1                      *
 *    pi = SUMMATORY( ------ [ ------  - ------ - ------ - ------]),  n >=0         *
 *                     16^n    8n + 1    8n + 4   8n + 5   8n + 6                   *
 *                                                                                  *
 * Formula quotients are coded as:                                                  *
 *              4                 2                 1                 1             *
 *   quot_a = ------,  quot_b = ------,  quot_c = ------,  quot_d = ------,         *
 *            8n + 1            8n + 4            8n + 5            8n + 6          *
 *                                                                                  *
 *              1                                                                   *
 *   quot_m = -----                                                                 *
 *             16^n                                                                 *
 *                                                                                  *
 ************************************************************************************
 * Bailey Borwein Plouffe formula dependencies:                                     *
 *                                                                                  *
 *                        1            1                                            *
 *           dep_m(n) = ----- = ---------------                                     *
 *                       16^n   dep_m(n-1) * 16                                     *
 *                                                                                  *
 ************************************************************************************/


void bbp_blocks_algorithm_gmp(mpf_t pi, int num_iterations, int num_threads){
    mpf_t quotient; 

    mpf_init_set_d(quotient, QUOTIENT);         // quotient = (1 / 16)   

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, block_size, block_start, block_end;
        mpf_t local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux;

        thread_id = omp_get_thread_num();
        block_size = (num_iterations + num_threads - 1) / num_threads;
        block_start = thread_id * block_size;
        block_end = block_start + block_size;
        if (block_end > num_iterations) block_end = num_iterations;

        mpf_init_set_ui(local_pi, 0);               // private thread pi
        mpf_init(dep_m);
        mpf_pow_ui(dep_m, quotient, block_start);    // m = (1/16)^n                  
        mpf_inits(quot_a, quot_b, quot_c, quot_d, aux, NULL);

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = block_start; i < block_end; i++){
                bbp_iteration_gmp(local_pi, i, dep_m, quot_a, quot_b, quot_c, quot_d, aux);
                // Update dependencies:  
                mpf_mul(dep_m, dep_m, quotient);
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpf_add(pi, pi, local_pi);

        //Clear thread memory
        mpf_clears(local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux, NULL);   
    }
        
    //Clear memory
    mpf_clear(quotient);
}
