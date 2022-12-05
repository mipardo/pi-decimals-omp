#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "chudnovsky_blocks_with_simplified_expression.h"


#define A 13591409
#define B 545140134
#define C 640320
#define D 426880
#define E 10005


/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Chudnovsky formula implementation                                                *
 * This version does not computes all the factorials (the expression is simplified) *
 * This version uses a "non-proportional" by blocks distribution                    *
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
 * This method provides an optimal distribution for each thread
 * based on the Chudnovsky iterations analysis.
 * It returns an array of three integers:
 *   distribution[0] -> block size
 *   distribution[1] -> block start
 *   distribution[2] -> block end 
 */
int * get_thread_distribution(int num_threads, int thread_id, int num_iterations){
    int * distribution, i, block_size, block_start, block_end, row, column; 
    FILE * ratios_file;
    float working_ratios[160][41];

    //Open the working_ratios file 
    ratios_file = fopen("resources/working_ratios.txt", "r");
    if(ratios_file == NULL){
        printf("working_ratios.txt not found \n");
        exit(-1);
    } 

    //Load the working_ratios matrix 
    row = 0;
    while (fscanf(ratios_file, "%f", &working_ratios[row][0]) == 1){
        for (column = 1; column < 41; column++){
            fscanf(ratios_file, "%f", &working_ratios[row][column]);
        }
        row++;
    }

    distribution = malloc(sizeof(int) * 3);
    if(num_threads == 1){
        distribution[0] = num_iterations;
        distribution[1] = 0;
        distribution[2] = num_iterations;
        return distribution; 
    }

    block_size = working_ratios[thread_id][num_threads / 4] * num_iterations / 100;
    block_start = 0;
    for(i = 0; i < thread_id; i ++){
        block_start += working_ratios[i][num_threads / 4] * num_iterations / 100;
    }
    block_end = block_start + block_size;
    if (thread_id == num_threads -1) block_end = num_iterations;
    
    distribution[0] = block_size;
    distribution[1] = block_start;
    distribution[2] = block_end;

    return distribution;
}



void chudnovsky_cheater_with_simplified_expression_algorithm_gmp(mpf_t pi, int num_iterations, int num_threads){
    mpf_t e, c;

    mpf_init_set_ui(e, E);
    mpf_init_set_ui(c, C);
    mpf_neg(c, c);
    mpf_pow_ui(c, c, 3);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {   
        int thread_id, i, block_size, block_start, block_end, factor_a, * distribution;
        mpf_t local_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux;

        thread_id = omp_get_thread_num();
        distribution = get_thread_distribution(num_threads, thread_id, num_iterations);
        block_size = distribution[0];
        block_start = distribution[1];
        block_end = distribution[2];
        
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