#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Print_Title.h"
#include "../GMP/Pi_Calculator_GMP.h"
#include "../MPFR/Pi_Calculator_MPFR.h"



int incorrect_params(char* exec_name){
    printf("  Number of params are not correct. Try with:\n");
    printf("    %s library algorithm precision num_threads \n", exec_name);
    printf("\n");
}

int main(int argc, char **argv){    

    print_title();
    printf("  Parallel OMP version! \n");
    printf("\n");

    //Check the number of parameters are correct
    if(argc != 5){
        incorrect_params(argv[0]);
        exit(-1);
    }

    //Take algorithm and precision from params
    char *library = argv[1];
    int algorithm = atoi(argv[2]);    
    int precision = atoi(argv[3]);
    int num_threads = atoi(argv[4]);

    if (strcmp(library, "GMP") == 0) {
        calculate_pi_gmp(algorithm, precision, num_threads);
    } 
    else if (strcmp(library, "MPFR") == 0) {
        calculate_pi_mpfr(algorithm, precision, num_threads);
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
