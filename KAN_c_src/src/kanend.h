#ifndef KANEND_H
#define KANEND_H

#include <stdlib.h>
#include "urysohn.h"
typedef struct KANAddend{
    double muInner;
    double muOuter;
    double targetMin;
    double targetMax;
    double lastInnerValue;
    Urysohn* u;
    Univariate* univariate;
} KANAddend;

KANAddend* KANAddend_create(double* xmin, double* xmax, 
    double targetMin, double targetMax,
    int inner, int outer, double muInner, double muOuter, int number_of_inputs) ;

void KANAddend_destroy(KANAddend* addend) ;

KANAddend* KANAddend_copy(const KANAddend* addend) ;

void KANAddend_UpdateUsingMemory(KANAddend* addend, double diff) ;

void KANAddend_UpdateUsingInput(KANAddend* addend, double* input, double diff) ;

double KANAddend_ComputeUsingInput(KANAddend* addend, double* input) ;

#endif