#ifndef univariate_h
#define univariate_h

#include <stdlib.h>

typedef struct Univariate{
    int points;
    double c;
    double xmin;
    double xmax;
    double deltax;
    double *y;
    int lastLeftIndex;
    double lastLeftOffset;
} Univariate;


double relu(double x) ;

double relu_derivative(double x) ;
void Univariate_Init(Univariate *uni, double xmin, double xmax, double ymin, double ymax, int points,int size) ;

void Univariate_Destroy(Univariate *uni) ;

void Univariate_Copy(Univariate *dest, const Univariate *src) ;
void Univariate_SetLimits(Univariate *uni) ;

void Univariate_SetRandomFunction(Univariate *uni, double ymin, double ymax) ;

void Univariate_FitDefinition(Univariate *uni, double x) ;

double Univariate_GetDerivative(Univariate *uni, double x) ;

void Univariate_UpdateUsingInput(Univariate *uni, double x, double delta, double mu) ;

void Univariate_UpdateUsingMemory(Univariate *uni, double delta, double mu) ;

double Univariate_GetFunctionUsingInput(Univariate *uni, double x) ;

#endif