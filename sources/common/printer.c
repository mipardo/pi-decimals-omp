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
    file = fopen("resources/pi_decimals_title.txt", "r");
    if(file == NULL){
        printf("pi_decimals_title.txt not found \n");
        exit(-1);
    } 

    char character;
    int i = 0;
    while((character = fgetc(file)) != EOF){
        printf("%c", character);
    }
    
    fclose(file);
}

void check_errors(int precision, int num_iterations, int num_threads){
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

void print_results(char *library, char *algorithm_tag, int precision, int num_iterations, int num_threads, int decimals_computed, double execution_time) {
    printf("  Library used: %s \n", library);
    printf("  Algorithm: %s \n", algorithm_tag);
    printf("  Precision used: %d \n", precision);
    printf("  Number of iterations: %d \n", num_iterations);
    printf("  Number of threads: %d \n", num_threads);
    if (decimals_computed >= precision) { printf("  Correct decimals: %d \n", decimals_computed); } 
    else { printf("  Something went wrong. The execution just achieved %d decimals \n", decimals_computed); }
    printf("  Execution time: %f seconds \n", execution_time);
    printf("\n");
}

void print_results_csv(char *library, char *algorithm_tag, int precision, int num_iterations, int num_threads, int decimals_computed, double execution_time) {
    //char execution_time_string[100];
    //replace_decimal_point_by_coma(execution_time, execution_time_string);
    printf("OMP;");
    printf("%s;", library);
    printf("%s;", algorithm_tag);
    printf("%d;", precision);
    printf("%d;", num_iterations);
    printf("%d;", num_threads);
    printf("%d;", decimals_computed);
    printf("%f;\n", execution_time);
}

