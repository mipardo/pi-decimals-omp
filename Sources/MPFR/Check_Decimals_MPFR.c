#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>


int check_decimals(mpfr_t pi){
    //Cast the number we want to check to string
    int bytes_of_pi =  ( (int) ceil((float) pi -> _mpfr_prec / (float) GMP_NUMB_BITS) ) * sizeof(mp_limb_t);
    char calculated_pi[bytes_of_pi * 8]; 
    mpfr_sprintf(calculated_pi, "%Re", pi);

    //Read the correct pi number from numeroPiCorrecto.txt file and compares the decimals to calculated pi
    FILE * file;
    file = fopen("Resources/numeroPiCorrecto.txt", "r");
    if(file == NULL){
        printf("numeroPiCorrecto.txt not found \n");
        exit(-1);
    } 

    char correct_pi_char;
    int i = 0;
    while((correct_pi_char = fgetc(file)) != EOF){
        if( (i >= (bytes_of_pi * 8)) || (correct_pi_char != calculated_pi[i])){
            break;
        }
        i++;
    }
    i = (i < 2) ? 0: i - 2;
    
    fclose(file);

    //mpfr_printf("Pi: %Re \n", pi);

    return i;
}