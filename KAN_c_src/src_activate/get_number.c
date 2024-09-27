// The following code is translated from C++ to C.
// Note: This translation assumes the presence of necessary header files for math functions and dynamic memory allocation.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Formula3{
    double* target;
    double** inputs;
    int _N;
    int nInputs;
} Formula3;

void Formula3_Destroy(Formula3* formula) {
    free(formula->target);
    for (int i = 0; i < formula->_N; ++i) {
        free(formula->inputs[i]);
    }
    free(formula->inputs);
}

double Formula3_Function(Formula3* formula, double* input) {
    // y = (1/pi)*(2+2*x3)*(1/3)*(atan(20*exp(x5)*(x1-0.5+x2/6))+pi/2) + (1/pi)*(2+2*x4)*(1/3)*(atan(20*exp(x5)*(x1-0.5-x2/6))+pi/2);
    double pi = 3.14159265359;

    double y = (1.0 / pi);
    y *= (2.0 + 2.0 * input[2]);
    y *= (1.0 / 3.0);
    y *= atan(20.0 * (input[0] - 0.5 + input[1] / 6.0) * exp(input[4])) + pi / 2.0;

    double z = (1.0 / pi);
    z *= (2.0 + 2.0 * input[3]);
    z *= (1.0 / 3.0);
    z *= atan(20.0 * (input[0] - 0.5 - input[1] / 6.0) * exp(input[4])) + pi / 2.0;

    return y + z;
}


void Formula3_GetInput(Formula3* formula, double* input) {
    for (int i = 0; i < formula->nInputs; ++i) {
        input[i] = (double)(rand() % 1000) / 1000.0;
    }
}

double Formula3_GetTarget(Formula3* formula, double* input) {
    return Formula3_Function(formula, input);
}

void Formula3_GenerateData(Formula3* formula, int N) {
    formula->_N = N;
    formula->target = (double*)malloc(N * sizeof(double));
    if (formula->target == NULL) {
        perror("Failed to allocate memory for target");
        exit(EXIT_FAILURE);
    }
    formula->inputs = (double**)malloc(N * sizeof(double*));
    if (formula->inputs == NULL) {
        perror("Failed to allocate memory for inputs");
        free(formula->target);
        exit(EXIT_FAILURE);
    }
    int counter = 0;
    while (counter<N) {
        formula->inputs[counter] = (double*)malloc(formula->nInputs * sizeof(double));
        if (formula->inputs[counter] == NULL) {
            perror("Failed to allocate memory for inputs[counter]");
            // free memory
            for (int j = 0; j < counter; ++j) {
                free(formula->inputs[j]);
            }
            }
        Formula3_GetInput(formula, formula->inputs[counter]);
        formula->target[counter] = Formula3_Function(formula, formula->inputs[counter]);
        if (++counter >= N) break;
    }
}