#ifndef URYSOHN_H
#define URYSOHN_H

#include <stdlib.h>
#include <stdio.h>


typedef struct Univariate Univariate; // Assuming Univariate is defined elsewhere

typedef struct Urysohn {
    int length;
    Univariate** univariateList;
} Urysohn;

Urysohn* Urysohn_Create(double* xmin, double* xmax, double targetMin, double targetMax,
                        int* layers, int len) ;
void Urysohn_Destroy(Urysohn* instance) ;
Urysohn* Urysohn_Copy(const Urysohn* uri) ;

void Urysohn_UpdateUsingInput(Urysohn* instance, double delta, double* inputs, double mu) ;

void Urysohn_UpdateUsingMemory(Urysohn* instance, double delta, double mu) ;
double Urysohn_GetValueUsingInput(Urysohn* instance, double* inputs) ;


#endif