#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>

#define QUOTIENT 0.0625

/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Last version of Bailey Borwein Plouffe formula implementation                    *
 * It allows to compute pi using multiple threads                                   *
 * It uses a cyclic distribution                                                    *
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
 *                        1                   1                                     *
 *           dep_m(n) = ----- = -------------------------------------               *
 *                       16^n   dep_m(n-num_threads) * 16^num_threads               *
 *                                                                                  *
 ************************************************************************************/

/*
 * An iteration of Bailey Borwein Plouffe formula
 */
void bbp_iteration_gmp(mpf_t pi, int n, mpf_t dep_m, mpf_t quot_a, mpf_t quot_b, mpf_t quot_c, mpf_t quot_d, mpf_t aux){
    mpf_set_ui(quot_a, 4);              // quot_a = ( 4 / (8n + 1))
    mpf_set_ui(quot_b, 2);              // quot_b = (-2 / (8n + 4))
    mpf_set_ui(quot_c, 1);              // quot_c = (-1 / (8n + 5))
    mpf_set_ui(quot_d, 1);              // quot_d = (-1 / (8n + 6))
    mpf_set_ui(aux, 0);                 // aux = a + b + c + d  

    int i = n << 3;                     // i = 8n
    mpf_div_ui(quot_a, quot_a, i | 1);  // 4 / (8n + 1)
    mpf_div_ui(quot_b, quot_b, i | 4);  // 2 / (8n + 4)
    mpf_div_ui(quot_c, quot_c, i | 5);  // 1 / (8n + 5)
    mpf_div_ui(quot_d, quot_d, i | 6);  // 1 / (8n + 6)

    // aux = (a - b - c - d)   
    mpf_sub(aux, quot_a, quot_b);
    mpf_sub(aux, aux, quot_c);
    mpf_sub(aux, aux, quot_d);

    // aux = m * aux 
    mpf_mul(aux, aux, dep_m);   
    
    mpf_add(pi, pi, aux);  
}


void bbp_cyclic_algorithm_gmp(mpf_t pi, int num_iterations, int num_threads){
    mpf_t jump, quotient; 

    mpf_init_set_d(quotient, QUOTIENT);         // quotient = (1 / 16)   
    mpf_init_set_ui(jump, 1);        
    mpf_pow_ui(jump, quotient, num_threads);    // jump = (1/16)^num_threads

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i;
        mpf_t local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux;

        thread_id = omp_get_thread_num();
        mpf_init_set_ui(local_pi, 0);               // private thread pi
        mpf_init(dep_m);
        mpf_pow_ui(dep_m, quotient, thread_id);    // m = (1/16)^n                  
        mpf_inits(quot_a, quot_b, quot_c, quot_d, aux, NULL);

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = thread_id; i < num_iterations; i+=num_threads){
                bbp_iteration_gmp(local_pi, i, dep_m, quot_a, quot_b, quot_c, quot_d, aux);
                // Update dependencies:  
                mpf_mul(dep_m, dep_m, jump);
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpf_add(pi, pi, local_pi);

        //Clear thread memory
        mpf_clears(local_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux, NULL);   
    }
        
    //Clear memory
    mpf_clears(jump, quotient, NULL);
}