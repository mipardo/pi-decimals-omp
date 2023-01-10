#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <time.h>
#include <stdbool.h>
#include "../common/printer.h"
#include "check_decimals.h"
#include "algorithms/bbp_blocks.h"
#include "algorithms/bellard_cyclic.h"
#include "algorithms/chudnovsky_simplified_expression_blocks.h"
#include "algorithms/chudnovsky_craig_wood_expression.h"


double gettimeofday();


void calculate_pi_mpfr(int algorithm, int precision, int num_threads, bool print_in_csv_format){
    double execution_time;
    struct timeval t1, t2;
    mpfr_t pi;
    int num_iterations, decimals_computed, precision_bits;
    char *algorithm_tag;

    precision_bits = 8 * precision;
    
    gettimeofday(&t1, NULL);

    //Set mpfr float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "MPFR-BBP-BLC";
        mpfr_bbp_blocks_algorithm(pi, num_iterations, num_threads, precision_bits);
        break;

    case 1:
        num_iterations = precision / 3;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "MPFR-BEL-CYC";
        mpfr_bellard_cyclic_algorithm(pi, num_iterations, num_threads, precision_bits);
        break;

    case 2:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "MPFR-CHD-SME-BLC";
        mpfr_chudnovsky_simplified_expression_blocks_algorithm(pi, num_iterations, num_threads, precision_bits);
        break;
    
    case 3:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "MPFR-CHD-CWE-SEQ";
        mpfr_chudnovsky_craig_wood_expression_algorithm(pi, num_iterations, num_threads, precision_bits);
        break;
    
    default:
        printf("  Algorithm number selected not available, try with another number. \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = mpfr_check_decimals(pi);
    if (print_in_csv_format) { print_results_csv("MPFR", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time); } 
    else { print_results("MPFR", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time); }
    mpfr_clear(pi);
}

