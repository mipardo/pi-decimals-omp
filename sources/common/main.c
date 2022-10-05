#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "printer.h"
#include "../gmp/pi_calculator.h"
#include "../mpfr/pi_calculator.h"



int incorrect_params(char* exec_name){
    printf("  Number of params are not correct. Try with:\n");
    printf("    %s library algorithm precision num_threads [-csv] \n", exec_name);
    printf("\n");
}

int main(int argc, char **argv){    

    //Check the number of parameters are correct
    bool print_in_csv_format; 

    if (argc == 6 && strcmp(argv[5], "-csv") == 0) {
        print_in_csv_format = true;    
    } 
    else if (argc == 5) {
        print_title();
        print_in_csv_format = false;
    }
    else {
        incorrect_params(argv[0]);
        exit(-1);
    }

    //Take algorithm and precision from params
    char *library = argv[1];
    int algorithm = atoi(argv[2]);    
    int precision = atoi(argv[3]);
    int num_threads = (atoi(argv[4]) <= 0) ? 1 : atoi(argv[4]);

    if (strcmp(library, "GMP") == 0) {
        calculate_pi_gmp(algorithm, precision, num_threads, print_in_csv_format);
    } 
    else if (strcmp(library, "MPFR") == 0) {
        calculate_pi_mpfr(algorithm, precision, num_threads, print_in_csv_format);
    } 
    else 
    {
        printf("  Library selected is not correct. Try with: \n");
        printf("      GMP -> GNU Multiple Precision Arithmetic Library \n");
        printf("      MPFR -> Multiple Precision Floating Point Reliable Library \n");
        printf("\n");
        exit(-1);
    }

    exit(0);
}
