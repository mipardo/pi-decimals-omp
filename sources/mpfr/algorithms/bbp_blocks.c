#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>


#define QUOTIENT 0.0625

/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Last version of Bailey Borwein Plouffe formula implementation                    *
 * It implements a single-threaded method and another that can use multiple threads *
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

/*
 * An iteration of Bailey Borwein Plouffe formula
 */
void bbp_iteration_mpfr(mpfr_t pi, int n, mpfr_t dep_m, mpfr_t quot_a, mpfr_t quot_b, mpfr_t quot_c, mpfr_t quot_d, mpfr_t aux){
    mpfr_set_ui(quot_a, 4, MPFR_RNDN);              // quot_a = ( 4 / (8n + 1))
    mpfr_set_ui(quot_b, 2, MPFR_RNDN);              // quot_b = (-2 / (8n + 4))
    mpfr_set_ui(quot_c, 1, MPFR_RNDN);              // quot_c = (-1 / (8n + 5))
    mpfr_set_ui(quot_d, 1, MPFR_RNDN);              // quot_d = (-1 / (8n + 6))
    mpfr_set_ui(aux, 0, MPFR_RNDN);                 // aux = a + b + c + d  

    int i = n << 3;                     // i = 8n
    mpfr_div_ui(quot_a, quot_a, i | 1, MPFR_RNDN);  // 4 / (8n + 1)
    mpfr_div_ui(quot_b, quot_b, i | 4, MPFR_RNDN);  // 2 / (8n + 4)
    mpfr_div_ui(quot_c, quot_c, i | 5, MPFR_RNDN);  // 1 / (8n + 5)
    mpfr_div_ui(quot_d, quot_d, i | 6, MPFR_RNDN);  // 1 / (8n + 6)

    // aux = (a - b - c - d)   
    mpfr_sub(aux, quot_a, quot_b, MPFR_RNDN);
    mpfr_sub(aux, aux, quot_c, MPFR_RNDN);
    mpfr_sub(aux, aux, quot_d, MPFR_RNDN);

    // aux = m * aux 
    mpfr_mul(aux, aux, dep_m, MPFR_RNDN);   
    
    mpfr_add(pi, pi, aux, MPFR_RNDN);  
}


void bbp_blocks_algorithm_mpfr(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
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
                bbp_iteration_mpfr(local_pi, i, dep_m, quot_a, quot_b, quot_c, quot_d, aux);
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