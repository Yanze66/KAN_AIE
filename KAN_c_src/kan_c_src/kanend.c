// This code is related to a custom mathematical or computational framework.

#include "kanend.h"
#include <stdlib.h>
#include "urysohn.h"
#include "univariat.h"
// typedef struct {
//     double muInner;
//     double muOuter;
//     double targetMin;
//     double targetMax;
//     double lastInnerValue;
//     Urysohn* u;
//     Univariate* univariate;
// } KANAddend;

KANAddend* KANAddend_create(double* xmin, double* xmax, 
    double targetMin, double targetMax,
    int inner, int outer, double muInner, double muOuter, int number_of_inputs) {

    KANAddend* addend = (KANAddend*)malloc(sizeof(KANAddend));
    addend->muInner = muInner;
    addend->muOuter = muOuter;
    addend->targetMin = targetMin;
    addend->targetMax = targetMax;
    addend->lastInnerValue = 0.0;

    int* interior_structure = (int*)malloc(number_of_inputs * sizeof(int));
    for (int i = 0; i < number_of_inputs; i++) {
        interior_structure[i] = inner;
    }
    
    addend->u = (Urysohn*)malloc(sizeof(Urysohn));
    addend->u = Urysohn_Create(xmin, xmax, targetMin, targetMax, interior_structure, number_of_inputs);
    addend->univariate = (Univariate*)malloc(sizeof(Univariate));
    Univariate_Init(addend->univariate,targetMin, targetMax, targetMin, targetMax, outer);
    
    free(interior_structure);
    return addend;
}

void KANAddend_destroy(KANAddend* addend) {
    if (addend) {
        Univariate_Destroy(addend->univariate);
        Urysohn_Destroy(addend->u);
        free(addend);
    }
}

KANAddend* KANAddend_copy(const KANAddend* addend) {
    KANAddend* newAddend = (KANAddend*)malloc(sizeof(KANAddend));
    newAddend->muInner = addend->muInner;
    newAddend->muOuter = addend->muOuter;
    newAddend->targetMin = addend->targetMin;
    newAddend->targetMax = addend->targetMax;
    newAddend->lastInnerValue = addend->lastInnerValue;
    newAddend->univariate = (Univariate*)malloc(sizeof(Univariate));
    Univariate_Copy(newAddend->univariate,addend->univariate);
    newAddend->u = (Urysohn*)malloc(sizeof(Urysohn));
    newAddend->u = Urysohn_Copy(addend->u);
    return newAddend;
}

void KANAddend_UpdateUsingMemory(KANAddend* addend, double diff) {
    double derivative = Univariate_GetDerivative(addend->univariate, addend->lastInnerValue);
    Urysohn_UpdateUsingMemory(addend->u, diff * derivative, addend->muInner);
    Univariate_UpdateUsingMemory(addend->univariate, diff, addend->muOuter);
}

void KANAddend_UpdateUsingInput(KANAddend* addend, double* input, double diff) {
    double value = Urysohn_GetValueUsingInput(addend->u, input);
    double derivative = Univariate_GetDerivative(addend->univariate, value);
    Urysohn_UpdateUsingInput(addend->u, diff * derivative, input, addend->muInner);
    Univariate_UpdateUsingInput(addend->univariate, value, diff, addend->muOuter);
}

double KANAddend_ComputeUsingInput(KANAddend* addend, double* input) {
    addend->lastInnerValue = Urysohn_GetValueUsingInput(addend->u, input);
    return Univariate_GetFunctionUsingInput(addend->univariate, addend->lastInnerValue);
}