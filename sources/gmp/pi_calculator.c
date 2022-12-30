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
#include "algorithms/chudnovsky_all_factorials_blocks.h"
#include "algorithms/chudnovsky_simplified_expression_blocks.h"
#include "algorithms/chudnovsky_simplified_expression_cheater.h"
#include "algorithms/chudnovsky_simplified_expression_snake_like.h"
#include "algorithms/chudnovsky_simplified_expression_integers_blocks.h"
#include "algorithms/chudnovsky_craig_wood_expression_blocks.h"


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
        gmp_bbp_cyclic_algorithm(pi, num_iterations, num_threads);
        break;

    case 1:
        num_iterations = precision * 0.84;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-BBP-BLC";
        gmp_bbp_blocks_algorithm(pi, num_iterations, num_threads);
        break;

    case 2:
        num_iterations = precision / 3;
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-BEL-CYC";
        gmp_bellard_cyclic_algorithm(pi, num_iterations, num_threads);
        break;

    case 3:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-CAF-BLC";
        gmp_chudnovsky_all_factorials_blocks_algorithm(pi, num_iterations, num_threads);
        break;

    case 4:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-SME-BLC";
        gmp_chudnovsky_simplified_expression_blocks_algorithm(pi, num_iterations, num_threads);
        break;

    case 5:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-SME-CHT";
        gmp_chudnovsky_simplified_expression_cheater_algorithm(pi, num_iterations, num_threads);
        break;

    case 6:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-SME-SNK";
        gmp_chudnovsky_simplified_expression_snake_like_algorithm(pi, num_iterations, num_threads);
        break;

    case 7:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-SME-INT-BLC";
        gmp_chudnovsky_simplified_expression_integers_blocks_algorithm(pi, num_iterations, num_threads);
        break;

    case 8:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations, num_threads);
        algorithm_tag = "GMP-CHD-CWE-BLC";
        gmp_chudnovsky_craig_wood_expression_blocks_algorithm(pi, num_iterations, num_threads);
        break;

    default:
        printf("  Algorithm number selected not availabe, try with another number. \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = gmp_check_decimals(pi);
    if (print_in_csv_format) {
        print_results_csv("GMP", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time);
    } else {
        print_results("GMP", algorithm_tag, precision, num_iterations, num_threads, decimals_computed, execution_time);
    }

    // gmp_printf("%.Ff \n", pi);
    mpf_clear(pi);

}