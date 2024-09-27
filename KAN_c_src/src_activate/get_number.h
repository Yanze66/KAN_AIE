#ifndef GET_NUMBER_H
#define GET_NUMBER_H
#include <math.h>
#include <stdlib.h>

typedef struct {
    double* target;
    double** inputs;
    int _N;
    int nInputs;
} Formula3;

void Formula3_Destroy(Formula3* formula) ;

double Formula3_Function(Formula3* formula, double* input) ;

void Formula3_GetInput(Formula3* formula, double* input) ;
double Formula3_GetTarget(Formula3* formula, double* input) ;
void Formula3_GenerateData(Formula3* formula, int N) ;
#endif