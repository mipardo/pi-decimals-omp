#include <stdio.h>
#include <stdlib.h>

void print_title(){
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
