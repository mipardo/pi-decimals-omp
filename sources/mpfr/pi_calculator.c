#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <time.h>
#include <stdbool.h>
#include "../common/printer.h"
#include "check_decimals.h"
#include "algorithms/bbp.h"
#include "algorithms/bellard.h"
#include "algorithms/chudnovsky.h"


double gettimeofday();

void check_errors_mpfr(int precision, int num_iterations, int num_threads){
    if (precision <= 0){
        printf("  Precision should be greater than cero. \n\n");
        exit(-1);
    } 
    if (num_iterations < num_threads){
        printf("  The number of iterations required for the computation is too small to be solved with %d threads. \n", num_threads);
        printf("  Try using a greater precision or lower threads number. \n\n");
        exit(-1);
    }
}

void calculate_pi_mpfr(int algorithm, int precision, int num_threads, bool print_in_csv_format){
    double execution_time;
    struct timeval t1, t2;
    mpfr_t pi;
    int num_iterations, decimals_computed, precision_bits;
    char *algorithm_type;

    precision_bits = 8 * precision;
    
    gettimeofday(&t1, NULL);

    //Set mpfr float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors_mpfr(precision, num_iterations, num_threads);
        algorithm_type = "BBP (Block distribution)";
        bbp_algorithm_mpfr(pi, num_iterations, num_threads, precision_bits);
        break;

    case 1:
        num_iterations = precision / 3;
        check_errors_mpfr(precision, num_iterations, num_threads);
        algorithm_type = "Bellard (Cyclic Distribution)";
        bellard_algorithm_mpfr(pi, num_iterations, num_threads, precision_bits);
        break;

    case 2:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors_mpfr(precision, num_iterations, num_threads);
        algorithm_type = "Chudnovsky (Block distribution and using the simplified mathematical expresion)";
        chudnovsky_algorithm_mpfr(pi, num_iterations, num_threads, precision_bits);
        break;
    
    default:
        printf("  Algorithm selected is not correct. Try with: \n");
        printf("      algorithm == 0 -> BBP (Block distribution) \n");
        printf("      algorithm == 1 -> Bellard (Cyclic Distribution) \n");
        printf("      algorithm == 2 -> Chudnovsky (Block distribution and using the simplified mathematical expresion) \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = check_decimals_mpfr(pi);
    if (print_in_csv_format) {
        print_results_csv("MPFR", algorithm_type, precision, num_iterations, num_threads, decimals_computed, execution_time);
    } 
    else {
        print_results("MPFR", algorithm_type, precision, num_iterations, num_threads, decimals_computed, execution_time);
    }
    mpfr_clear(pi);
}

