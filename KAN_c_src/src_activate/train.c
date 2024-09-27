// Framework: C, Standard Library
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "kanend.h"
#include "get_number.h"
#include <time.h>

void FindMinMax(double* xmin, double* xmax,
    double* targetMin, double* targetMax,
    double** matrix, double* target, int nRows, int nCols) {

    for (int columnIndex = 0; columnIndex < nCols; ++columnIndex) {
        xmin[columnIndex] = DBL_MAX;
        xmax[columnIndex] = -DBL_MIN;
    }

    for (int rowIndex = 0; rowIndex < nRows; ++rowIndex) {
        for (int columnIndex = 0; columnIndex < nCols; ++columnIndex) {
            if (matrix[rowIndex][columnIndex] < xmin[columnIndex]) {
                xmin[columnIndex] = matrix[rowIndex][columnIndex];
            }
            if (matrix[rowIndex][columnIndex] > xmax[columnIndex]) {
                xmax[columnIndex] = matrix[rowIndex][columnIndex];
            }
        }
    }

    *targetMin = DBL_MAX;
    *targetMax = -DBL_MIN;
    

    for (int rowIndex = 0; rowIndex < nRows; ++rowIndex) {
        if (target[rowIndex] < *targetMin) {
            *targetMin = target[rowIndex];
        }
        if (target[rowIndex] > *targetMax) {
            *targetMax = target[rowIndex];
        }
    }
}


void Training(double** inputs, double* target, 
    KANAddend** addends, int nRecords, int nEpochs, int nModels,
    int marginStart, int marginEnd, double sensitivity) {

    double* residualError = (double*)malloc(nRecords * sizeof(double));
    for (int i = 0; i < nRecords; ++i) {
        residualError[i] = 0.0;
    }
    int end = nEpochs - marginEnd;
    for (int epoch = 0; epoch < nEpochs; ++epoch) {
        double error2 = 0.0;
        int cnt = 0;
        for (int i = 0; i < nRecords; ++i) {
            if (epoch >= marginStart && epoch < end && residualError[i] < sensitivity) continue;
            double residual = target[i];
            for (int j = 0; j < nModels; ++j) {
                residual -= KANAddend_ComputeUsingInput(addends[j],inputs[i]);
            }
            for (int j = 0; j < nModels; ++j) {
                KANAddend_UpdateUsingInput(addends[j],inputs[i],residual);
            }
            error2 += residual * residual;
            residualError[i] = fabs(residual);
            ++cnt;
        }
        if (0 == cnt) error2 = 0.0;
        else {
            error2 /= cnt;
            error2 = sqrt(error2);
        }
        printf("Training step %d, current RMSE %4.4f\n", epoch, error2);
    }
    free(residualError);
}

int main(int argc, char **argv) {
    int nRecords = 10000;
    // int nModels = 11;
    int nModels = 30; //number of splines in output layer
    int nEpochs = 36;
    // int nEpochs = 100;
    int marginStart = 6;
    int marginEnd = 6;
    double sensitivity = 0.06;
    int innerPoints = 6; //control points of spline in hidden layer
    int outerPoints = 12;   //control points of spline in output layer
    double muInner = 0.01;
    double muOuter = 0.01;

    // Assuming Formula3 is defined and has appropriate methods
    Formula3 *formula = (Formula3 *)malloc(sizeof(Formula3));

    formula->nInputs = 5; //input dimension
    Formula3_GenerateData(formula, nRecords);

    double *xmin = (double *)malloc(formula->nInputs * sizeof(double));
    double *xmax = (double *)malloc(formula->nInputs * sizeof(double));
    double targetMin = 0.0; // Placeholder for targetMin
    double targetMax = 0.0; // Placeholder for targetMax

    FindMinMax(xmin, xmax, &targetMin, &targetMax, formula->inputs, formula->target, formula->_N, formula->nInputs);

    double zmin = targetMin / nModels;
    double zmax = targetMax / nModels;
    KANAddend **addends = (KANAddend **)malloc(nModels * sizeof(KANAddend *));
    for (int i = 0; i < nModels; ++i) {
        addends[i] = (KANAddend *)malloc(sizeof(KANAddend));
        addends[i] = KANAddend_create( xmin, xmax, zmin, zmax, innerPoints, outerPoints, muInner, muOuter, formula->nInputs);
    }

    clock_t start_encoding = clock();
    srand((unsigned int)time(NULL));
    

    Training( formula->inputs, formula->target, addends, nRecords, nEpochs, nModels, marginStart, marginEnd, sensitivity );
    
    clock_t end_encoding = clock();
    printf("\nTime for training %2.3f sec.\n", (double)(end_encoding - start_encoding) / CLOCKS_PER_SEC);
   
    printf("\nnumber of Parameter is %d\n", nModels * (innerPoints + 1) * 6); //univariate * 6
    //////// Object copy test //////
    KANAddend **addendsCopy = (KANAddend **)malloc(nModels * sizeof(KANAddend *));
    for (int i = 0; i < nModels; ++i) {
         addendsCopy[i] = (KANAddend *)malloc(sizeof(KANAddend));
         if (addendsCopy[i] == NULL) {
             perror("Failed to allocate memory for addendsCopy[i]");
             exit(EXIT_FAILURE);
         }
         addendsCopy[i] = KANAddend_copy( addends[i]);
    }

    double error = 0.0;
    double error3 = 0.0;
    int NTests = 1000;
    double *test_input = (double *)malloc(formula->nInputs * sizeof(double));
    for (int i = 0; i < NTests; ++i) {
        Formula3_GetInput(formula, test_input);
        double test_target = Formula3_GetTarget(formula, test_input);

        double model1 = 0.0;
        for (int j = 0; j < nModels; ++j) {
            model1 += KANAddend_ComputeUsingInput(addends[j], test_input);
        }

        double model2 = 0.0;
        for (int j = 0; j < nModels; ++j) {
            model2 += KANAddend_ComputeUsingInput(addendsCopy[j], test_input);
        }

        error += (test_target - model1) * (test_target - model1);
        error3 += (test_target - model2) * (test_target - model2);
    }
    error /= NTests;
    error = sqrt(error);
    error /= (targetMax - targetMin);
    error3 /= NTests;
    error3 = sqrt(error3);
    error3 /= (targetMax - targetMin);
    printf("\nRelative RMSE for unseen data %f, RMSE for copy object %f\n", error, error3);

    // Free allocated memory
    free(xmin);
    free(xmax);
    for (int i = 0; i < nModels; ++i) {
        free(addends[i]);
        free(addendsCopy[i]);
    }
    free(addends);
    free(addendsCopy);
    free(formula);
    free(test_input);

    return 0;
}