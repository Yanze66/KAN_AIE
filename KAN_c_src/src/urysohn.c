// Framework: C++ to C translation of a class-based implementation
#include <stdlib.h>
#include <stdio.h>
#include "univariat.h" // Assuming Univariate is defined elsewhere
#include "urysohn.h"

typedef struct Univariate Univariate; // Assuming Univariate is defined elsewhere


Urysohn* Urysohn_Create(double* xmin, double* xmax, double targetMin, double targetMax,
                        int* layers, int len) {
    Urysohn* instance = (Urysohn*)malloc(sizeof(Urysohn));
    if (instance == NULL) {
        perror("Failed to allocate memory for Urysohn");
        return NULL;
    }

    instance->length = len;
    double ymin = targetMin / instance->length;
    double ymax = targetMax / instance->length;

    instance->univariateList = (Univariate**)malloc(sizeof(Univariate*) * instance->length);
    if (instance->univariateList == NULL) {
        perror("Failed to allocate memory for univariateList");
        free(instance);
        return NULL;
    }

    for (int i = 0; i < instance->length; ++i) {
          instance->univariateList[i] = (Univariate*)malloc(sizeof(Univariate)); // 为每个 Univariate 分配内存
          if (instance->univariateList[i] == NULL) {
            perror("Memory allocation for uni failed");
            exit(EXIT_FAILURE);
            }
          Univariate_Init(instance->univariateList[i],xmin[i], xmax[i], ymin, ymax, layers[i]);
    }
    
    return instance;
}

void Urysohn_Destroy(Urysohn* instance) {
    for (int i = 0; i < instance->length; ++i) {
        Univariate_Destroy(instance->univariateList[i]); // Assuming Univariate_Destroy is defined
    }
    free(instance->univariateList);
    free(instance);
}

Urysohn* Urysohn_Copy(const Urysohn* uri) {
    Urysohn* instance = (Urysohn*)malloc(sizeof(Urysohn));
    instance->length = uri->length;
    instance->univariateList = (Univariate**)malloc(sizeof(Univariate*) * instance->length);
    
    for (int i = 0; i < instance->length; ++i) {
        instance->univariateList[i] = (Univariate*)malloc(sizeof(Univariate)); // 分配内存
        Univariate_Copy(instance->univariateList[i], uri->univariateList[i]); // Assuming Univariate_Copy is defined
    }
    
    return instance;
}

void Urysohn_UpdateUsingInput(Urysohn* instance, double delta, double* inputs, double mu) {
    for (int i = 0; i < instance->length; ++i) {
        Univariate_UpdateUsingInput(instance->univariateList[i], inputs[i], delta, mu); // Assuming Univariate_UpdateUsingInput is defined
    }
}

void Urysohn_UpdateUsingMemory(Urysohn* instance, double delta, double mu) {
    for (int i = 0; i < instance->length; ++i) {
        Univariate_UpdateUsingMemory(instance->univariateList[i], delta, mu); // Assuming Univariate_UpdateUsingMemory is defined
    }
}

double Urysohn_GetValueUsingInput(Urysohn* instance, double* inputs) {
    double f = 0.0;
    // printf("length: %d\n", instance->length);
    for (int i = 0; i < instance->length; ++i) {
        if (instance->univariateList[i] == NULL) {
            perror("Urysohn_GetValueUsingInput:Memory allocation for uni failed");
            exit(EXIT_FAILURE);
            }
        f += Univariate_GetFunctionUsingInput(instance->univariateList[i], inputs[i]); // Assuming Univariate_GetFunctionUsingInput is defined
    }
    return f;
}