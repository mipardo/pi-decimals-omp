#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void replace_decimal_point_by_coma(double number, char *result){
    int i;
    sprintf(result, "%f", number);
    for (i = 0; i < strlen(result); i++) {
        if (result[i] == '.'){
            result[i] = ',';
            break;
        }
    }
}

void print_title() {
    FILE * file;
    file = fopen("Resources/piDecimalsTitle.txt", "r");
    if(file == NULL){
        printf("piDecimalsTitle.txt not found \n");
        exit(-1);
    } 

    char character;
    int i = 0;
    while((character = fgetc(file)) != EOF){
        printf("%c", character);
    }
    
    fclose(file);
}

void print_results(char *library, char *algorithm_type, int precision, int num_iterations, int num_threads, int decimals_computed, double execution_time) {
    printf("  Library used: %s \n", library);
    printf("  Algorithm: %s \n", algorithm_type);
    printf("  Precision used: %d \n", precision);
    printf("  Number of iterations: %d \n", num_iterations);
    printf("  Number of threads: %d \n", num_threads);
    printf("  Correct decimals: %d \n", decimals_computed);
    printf("  Execution time: %f seconds \n", execution_time);
    printf("\n");
}

void print_results_csv(char *library, char *algorithm_type, int precision, int num_iterations, int num_threads, int decimals_computed, double execution_time) {
    char execution_time_string[100];
    replace_decimal_point_by_coma(execution_time, execution_time_string);
    //printf("library;algorithm_type;precision;num_iterations;num_threads;decimals_computed;execution_time;\n");
    printf("%s;%s;%d;%d;%d;%d;%s;\n", library, algorithm_type, precision, num_iterations, num_threads, decimals_computed, execution_time_string);
}

