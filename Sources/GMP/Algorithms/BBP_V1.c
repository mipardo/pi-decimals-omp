#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>


#define QUOTIENT 0.0625

/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * First version of Bailey Borwein Plouffe formula implementation                   *
 * It allows to compute pi using multiple threads                                   *
 * This is and old and slow implementation, therefore it is not used (DEPRECATED)   *
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
 *   quot_m = ------                                                                *
 *             16^n                                                                 *
 *                                                                                  *
 ************************************************************************************

/*
 * An iteration of Bailey Borwein Plouffe formula
 */
void bbp_iteration_v1_gmp(mpf_t pi, int n, mpf_t quotient){
    mpf_t quot_a, quot_b, quot_c, quot_d, quot_m, aux;

    mpf_init_set_ui(quot_a, 4);         // quot_a = (  4 / 8n + 1))
    mpf_init_set_ui(quot_b, 2);         // quot_b = ( -2 / 8n + 4))
    mpf_init_set_ui(quot_c, 1);         // quot_c = ( -1 / 8n + 5))
    mpf_init_set_ui(quot_d, 1);         // quot_d = ( -1 / 8n + 6))
    mpf_init_set_ui(quot_m, 0);         // quot_m = (1/16)^n  
    mpf_init(aux);                      // aux = a + b + c + d  

    int i = n * 8;                 
    mpf_div_ui(quot_a, quot_a, i + 1);  // 4 / (8n + 1)
    mpf_div_ui(quot_b, quot_b, i + 4);  // 2 / (8n + 4)
    mpf_div_ui(quot_c, quot_c, i + 5);  // 1 / (8n + 5)
    mpf_div_ui(quot_d, quot_d, i + 6);  // 1 / (8n + 6)

    // aux = (a - b - c - d)   
    mpf_sub(aux, quot_a, quot_b);
    mpf_sub(aux, aux, quot_c);
    mpf_sub(aux, aux, quot_d);

    // aux = m * aux
    mpf_pow_ui(quot_m, quotient, n);    // (1/16)^n
    mpf_mul(aux, aux, quot_m);          // quot_m * aux
    
    mpf_add(pi, pi, aux); 

    mpf_clears(quot_a, quot_b, quot_c, quot_d, aux, quot_m, NULL);
}


/*
 * Parallel Pi number calculation using the BBP algorithm
 * The number of iterations is divided cyclically, 
 * so each thread calculates a part of Pi.  
 */
void bbp_algorithm_v1_gmp(mpf_t pi, int num_iterations, int num_threads){
    int thread_id, i;
    mpf_t quotient; 

    mpf_init_set_d(quotient, QUOTIENT); // quotient = (1 / 16)   

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel private(thread_id, i)
    {
        mpf_t local_pi;

        thread_id = omp_get_thread_num();
        mpf_init_set_ui(local_pi, 0);   // private thread pi
        
        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = thread_id; i < num_iterations; i+=num_threads){
                bbp_iteration_v1_gmp(local_pi, i, quotient);    
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpf_add(pi, pi, local_pi);

        //Clear thread memory
        mpf_clear(local_pi);   
    }
        
    //Clear memory
    mpf_clear(quotient);
}
