#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "Check_Decimals_GMP.h"
#include "Algorithms/BBP_Blocks.h"
#include "Algorithms/BBP_Cyclic.h"
#include "Algorithms/Bellard.h"
#include "Algorithms/Chudnovsky_V1.h"
#include "Algorithms/Chudnovsky_V2.h"
#include "Algorithms/Chudnovsky.h"




double gettimeofday();

void check_errors_gmp(int precision, int num_iterations, int num_threads, int algorithm){
    if (precision <= 0){
        printf("  Precision should be greater than cero. \n\n");
        exit(-1);
    } 
    if (num_iterations < num_threads){
        printf("  The number of iterations required for the computation is too small to be solved with %d threads. \n", num_threads);
        printf("  Try using a greater precision or lower threads number. \n\n");
        exit(-1);
    }
    if (algorithm == 5){ 
        // Last version of Chudnovksy is more efficient when threads and procs are 2 or multiples of four
        if (num_threads > 2 && num_threads % 4 != 0){
            printf("  The last version of Chudnovksy is not eficient with %d threads. \n", num_threads);
            printf("  Try using 2 threads or multiples of four (4, 8, 12, 16, ..) \n\n");
            exit(-1);
        } 
    }
}

void print_running_properties_gmp(int precision, int num_iterations, int num_threads){
    printf("  Library used: GMP \n");
    printf("  Precision used: %d \n", precision);
    printf("  Number of iterations: %d \n", num_iterations);
    printf("  Number of threads: %d\n", num_threads);
}

void calculate_pi_gmp(int algorithm, int precision, int num_threads){
    double execution_time;
    struct timeval t1, t2;
    mpf_t pi;
    int num_iterations, decimals_computed;

    gettimeofday(&t1, NULL);

    //Set gmp float precision (in bits) and init pi
    mpf_set_default_prec(precision * 8); 
    mpf_init_set_ui(pi, 0); 
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: BBP (Cyclic distribution) \n");
        print_running_properties_gmp(precision, num_iterations, num_threads);
        bbp_algorithm_cyclic_gmp(pi, num_iterations, num_threads);
        break;

    case 1:
        num_iterations = precision * 0.84;
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: BBP (Block distribution)\n");
        print_running_properties_gmp(precision, num_iterations, num_threads);      
        bbp_algorithm_blocks_gmp(pi, num_iterations, num_threads);
        break;

    case 2:
        num_iterations = precision / 3;
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: Bellard (Cyclic Distribution) \n");
        print_running_properties_gmp(precision, num_iterations, num_threads);
        bellard_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 3:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: Chudnovsky (Block Distribution and computing all factorials) \n");
        print_running_properties_gmp(precision, num_iterations, num_threads);
        chudnovsky_algorithm_v1_gmp(pi, num_iterations, num_threads);
        break;

    case 4:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: Chudnovsky (Block distribution and using the simplified mathematical expresion) \n");
        print_running_properties_gmp(precision, num_iterations, num_threads);
        chudnovsky_algorithm_v2_gmp(pi, num_iterations, num_threads);
        break;

    case 5:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors_gmp(precision, num_iterations, num_threads, algorithm);
        printf("  Algorithm: Chudnovsky (Non-proportional block distribution and using the simplified mathematical expresion) \n");
        print_running_properties_gmp(precision, num_iterations, num_threads);
        chudnovsky_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    default:
        printf("  Algorithm selected is not correct. Try with: \n");
        printf("      algorithm == 0 -> BBP (Cyclic distribution) \n");
        printf("      algorithm == 1 -> BBP (Blocks distribution) \n");
        printf("      algorithm == 2 -> Bellard (First version) \n");
        printf("      algorithm == 3 -> Chudnovsky (Block Distribution and computing all factorials)\n");
        printf("      algorithm == 4 -> Chudnovsky (Block distribution and using the simplified mathematical expresion) \n");
        printf("      algorithm == 5 -> Chudnovsky (Non-proportional block distribution and using the simplified mathematical expresion) \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = check_decimals_gmp(pi);
    mpf_clear(pi);
    printf("  Match the first %d decimals. \n", decimals_computed);
    printf("  Execution time: %f seconds. \n", execution_time);
    printf("\n");
}