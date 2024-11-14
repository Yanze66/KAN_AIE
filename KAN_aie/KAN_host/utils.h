#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double identity(double x);
double square(double x);
double cube(double x);
double power4(double x);
double power5(double x);
double reciprocal(double x);
double reciprocal_square(double x);
double reciprocal_cube(double x);
double reciprocal_x4(double x);
double reciprocal_x5(double x);
double sqrt_func(double x);
double cube_root(double x);
double reciprocal_sqrt(double x);
double exp_func(double x);
double log_func(double x);
double abs_func(double x);
double sin_func(double x);
double cos_func(double x);
double tanh_func(double x);
double sgn_func(double x);
double arctanh_func(double x);
double arccos_func(double x);
double arcsin_func(double x);
double tan_func(double x);
double arctan_func(double x);
double zero_func(double x);
double gauss_func(double x);

// Derivative of y = x, i.e., identity function
double identity_derivative(double x) ;
// Derivative of y = x^2
double square_derivative(double x) ;
// Derivative of y = x^3
double cube_derivative(double x) ;

// Derivative of y = x^4
double power4_derivative(double x) ;
// Derivative of y = x^5
double power5_derivative(double x) ;

// Derivative of y = 1/x
double reciprocal_derivative(double x) ;
// Derivative of y = 1/x^2
double reciprocal_square_derivative(double x) ;

// Derivative of y = 1/x^3
double reciprocal_cube_derivative(double x) ;
// Derivative of y = 1/x^4
double reciprocal_x4_derivative(double x) ;
// Derivative of y = 1/x^5
double reciprocal_x5_derivative(double x) ;
// Derivative of y = sqrt(x)
double sqrt_func_derivative(double x) ;

// Derivative of y = sqrt(x)^3
double cube_root_derivative(double x) ;

// Derivative of y = 1/sqrt(x)
double reciprocal_sqrt_derivative(double x) ;

// Derivative of y = exp(x)
double exp_func_derivative(double x) ;

// Derivative of y = log(x)
double log_func_derivative(double x) ;
// Derivative of y = abs(x)
double abs_func_derivative(double x) ;

// Derivative of y = sin(x)
double sin_func_derivative(double x) ;
// Derivative of y = cos(x)
double cos_func_derivative(double x) ;
// Derivative of y = tan(x)
double tan_func_derivative(double x) ;
// Derivative of y = tanh(x)
double tanh_func_derivative(double x) ;
// Derivative of y = sgn(x)
double sgn_func_derivative(double x) ;
// Derivative of y = arcsin(x)
double arcsin_func_derivative(double x) ;
// Derivative of y = arccos(x)
double arccos_func_derivative(double x) ;
// Derivative of y = arctan(x)
double arctan_func_derivative(double x) ;
// Derivative of y = arctanh(x)
double arctanh_func_derivative(double x) ;

// Derivative of y = 0
double zero_func_derivative(double x) ;
// Derivative of y = exp(-x^2), often called the Gaussian function
double gauss_func_derivative(double x) ;

void generate_data(double *inputs, double *outputs, int num_samples) ;
double calculate_rmse(const double *actual, const double *predicted, int size) ;
void standardize_data(double *data, int num_samples) ;
double normal_random(double sigma);

#endif