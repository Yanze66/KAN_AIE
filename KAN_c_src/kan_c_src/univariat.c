// Translated C code for Univariate class functionality
#include <stdio.h>
#include <stdlib.h>
#include "univariat.h"
// typedef struct {
//     int points;
//     double xmin;
//     double xmax;
//     double deltax;
//     double *y;
//     int lastLeftIndex;
//     double lastLeftOffset;
// } Univariate;

void Univariate_Init(Univariate *uni, double xmin, double xmax, double ymin, double ymax, int points) {
    if (uni == NULL) {
        printf("uni pointer is NULL");
        return; // 或者其他错误处理方式
    }
    uni->points = points;
    uni->xmin = xmin;
    uni->xmax = xmax;
    Univariate_SetLimits(uni);
    Univariate_SetRandomFunction(uni, ymin, ymax);
}

void Univariate_Destroy(Univariate *uni) {
    free(uni->y);
}

void Univariate_Copy(Univariate *dest, const Univariate *src) {
    dest->points = src->points;
    dest->xmin = src->xmin;
    dest->xmax = src->xmax;
    dest->deltax = (src->xmax - src->xmin) / (src->points - 1);
    dest->y = (double *)malloc(src->points * sizeof(double));
    for (int i = 0; i < src->points; i++) {
        dest->y[i] = src->y[i];
    }
    dest->lastLeftIndex = src->lastLeftIndex;
    dest->lastLeftOffset = src->lastLeftOffset;
}

void Univariate_SetLimits(Univariate *uni) {
    double range = uni->xmax - uni->xmin;
    uni->xmin -= 0.01 * range;
    uni->xmax += 0.01 * range;
    uni->deltax = (uni->xmax - uni->xmin) / (uni->points - 1);
}

void Univariate_SetRandomFunction(Univariate *uni, double ymin, double ymax) {
    uni->y = (double *)malloc(uni->points * sizeof(double));
    for (int i = 0; i < uni->points; ++i) {
        uni->y[i] = rand() % 100;
    }
    double min = uni->y[0];
    double max = uni->y[0];
    for (int i = 0; i < uni->points; ++i) {
        if (uni->y[i] < min) min = uni->y[i];
        if (uni->y[i] > max) max = uni->y[i];
    }
    if (min == max) max = min + 1.0;
    for (int i = 0; i < uni->points; ++i) {
        uni->y[i] = (uni->y[i] - min) / (max - min) * (ymax - ymin) + ymin;
    }
}

void Univariate_FitDefinition(Univariate *uni, double x) {
    if (x < uni->xmin) {
        uni->xmin = x;
        Univariate_SetLimits(uni);
    }
    if (x > uni->xmax) {
        uni->xmax = x;
        Univariate_SetLimits(uni);
    }
}

double Univariate_GetDerivative(Univariate *uni, double x) {
    if (uni->deltax == 0) {
        fprintf(stderr, "Invalid deltax, cannot divide by zero\n");
        return 1;
    }
    int low = (int)((x - uni->xmin) / uni->deltax);
    return (uni->y[low + 1] - uni->y[low]) / uni->deltax;
}

void Univariate_UpdateUsingInput(Univariate *uni, double x, double delta, double mu) {
    if (uni == NULL) {
        perror("Univariate_UpdateUsingInput :Univariate pointer is NULL");
        return ;
    }
    Univariate_FitDefinition(uni, x);
    delta *= mu;
    double offset = (x - uni->xmin) / uni->deltax;
    int left = (int)(offset);
    double leftx = offset - left;
    uni->y[left + 1] += delta * leftx;
    uni->y[left] += delta * (1.0 - leftx);
}

void Univariate_UpdateUsingMemory(Univariate *uni, double delta, double mu) {
    delta *= mu;
    uni->y[uni->lastLeftIndex + 1] += delta * uni->lastLeftOffset;
    uni->y[uni->lastLeftIndex] += delta * (1.0 - uni->lastLeftOffset);
}

double Univariate_GetFunctionUsingInput(Univariate *uni, double x) {
    if (uni == NULL) {
        perror("Univariate_GetFunctionUsingInput:Univariate pointer is NULL");
        return 1;
    }
    Univariate_FitDefinition(uni, x);
    double offset = (x - uni->xmin) / uni->deltax;
    int leftIndex = (int)(offset);
    double leftOffset = offset - leftIndex;
    uni->lastLeftIndex = leftIndex;
    uni->lastLeftOffset = leftOffset;
    return uni->y[leftIndex] + (uni->y[leftIndex + 1] - uni->y[leftIndex]) * leftOffset;
}