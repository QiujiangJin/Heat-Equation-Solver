/*********************************************************************
 * This is the main function of this project. It creates a solver of 
 * ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

    double res = 0.0;

    FILE *f_1 = fopen(argv[1], "r");
    if(f_1 == NULL)  {
        printf("Error opening file!\n");
        exit(1);
    }
    FILE *f_2 = fopen(argv[2], "r");
    if(f_2 == NULL)  {
        printf("Error opening file!\n");
        exit(1);
    }

    int line = 0;
    int ch = 0;

    while(!feof(f_1)) {
  	ch = fgetc(f_1);
  	if(ch == '\n'){
    	    line++;
  	}
    }

    double n_1 = 0.0;
    double n_2 = 0.0;

    for(int i = 0; i < line; i++) {
	fscanf(f_1, "%lf", &n_1);
	fscanf(f_2, "%lf", &n_2);
	res += (n_1 - n_2)*(n_1 - n_2);
    } 

    fclose(f_1);
    fclose(f_2);

    printf("%f", sqrt(res));

}
