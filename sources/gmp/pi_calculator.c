#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <stdbool.h>
#include "../common/printer.h"
#include "check_decimals.h"
#include "algorithms/bbp_blocks.h"
#include "algorithms/bbp_cyclic.h"
#include "algorithms/bellard_cyclic.h"
#include "algorithms/chudnovsky_blocks_with_all_factorials.h"
#include "algorithms/chudnovsky_blocks_with_simplified_expression.h"
#include "algorithms/chudnovsky_cheater_with_simplified_expression.h"
#include "algorithms/chudnovsky_snake_like_with_simplified_expression.h"


double gettimeofday();


void calculate_pi_gmp(int algorithm, int precision, int num_threads, bool print_in_csv_format){
    double execution_time;
    struct timeval t1, t2;
    mpf_t pi;
    int num_iterations, decimals_computed;
    char *algorithm_tag;

    gettimeofday(&t1, NULL);

    //Set gmp float precision (in bits) and init pi
    mpf_set_default_prec(precision * 8); 
    mpf_init_set_ui(pi, 0); 
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-BBP-CYC";
        bbp_cyclic_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 1:
        num_iterations = precision * 0.84;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-BBP-BLC";
        bbp_blocks_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 2:
        num_iterations = precision / 3;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-BEL-CYC";
        bellard_cyclic_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 3:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-BLC-CAF";
        chudnovsky_blocks_with_all_factorials_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 4:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-BLC-SME";
        chudnovsky_blocks_with_simplified_expression_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 5:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-CHT-SME";
        chudnovsky_cheater_with_simplified_expression_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    case 6:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-SNK-SME";
        chudnovsky_snake_like_with_simplified_expression_algorithm_gmp(pi, num_iterations, num_threads);
        break;

    default:
        printf("  Algorithm number selected not availabe, try with another number. \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = check_decimals_gmp(pi);
    if (print_in_csv_format) {
        print_results_csv("GMP", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time);
    } else {
        print_results("GMP", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time);
    }
    mpf_clear(pi);

}