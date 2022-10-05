#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>


int check_decimals_gmp(mpf_t pi){
    //Cast the number we want to check to string
    int bytes_of_pi = ((pi -> _mp_prec + 1) * sizeof(mp_limb_t)) * sizeof(int); 
    char calculated_pi[bytes_of_pi]; 
    gmp_sprintf(calculated_pi, "%.Ff", pi);

    //Read the correct pi number from numeroPiCorrecto.txt file and compares the decimals to calculated pi
    FILE * file;
    file = fopen("resources/numero_pi_correcto.txt", "r");
    if(file == NULL){
        printf("numero_pi_correcto.txt not found \n");
        exit(-1);
    } 

    char correct_pi_char;
    int i = 0;
    while((correct_pi_char = fgetc(file)) != EOF){
        if( (i >= bytes_of_pi) || (correct_pi_char != calculated_pi[i])){
            break;
        }
        i++;
    }
    i = (i < 2) ? 0: i - 2;
    
    fclose(file);

    return i;
}
