#ifndef SPLINE2FUNLOOP_H
#define SPLINE2FUNLOOP_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"

double calculate_bspline_derivative(double *grid, int k, int num_of_basis, double x) ;



//  double calculate_bspline(double *grid, int k, int num_of_basis, double x) ;
double calculate_bspline(double *grid, int num_of_basis, double x) ;

//  define the function to calculate the B-spline basis function,x为输入，grid为节点，num为节点数，k为阶数，result为输出
void B_batch(double x, double *grid, int num, int k, double *result) ;

double relu(double x) ;

// define the structure of a single neuron
typedef struct {
    int layer; // layer index
     
    int in_index; // input index
    int out_index; // output index

    double *output_history;  // 存储最近20次计算的输出
    int history_index;          // 当前索引
    


    double wb;   // weight of the bias
    int num; // spline nodes
    int k;   // spline order
    int grid_size; // grid size
    double *coef; // Bspline coefficients
    double *grid; // Bspline grid
    int active;   // neuron activation status

    double in_value; // neuron input value
    double out_value; // neuron output value

    int fun_fixed; // function fixed status
    double (*fx)(double);  // 函数指针用于动态确定输出
    double (*fx_derivative)(double);  // 函数指针用于动态确定输出
    double best_a; // scale factor
    double best_b; // shift factor
    double best_c; // bias factor
    double best_d; // slope factor

    double pre_error; // pre_error of the neuron
} Neural;

Neural* init_neural(int l, int in_index, int out_index, int num, int k) ;

double forward_neural(Neural *neural, double x) ;

//void grid_extend(double x, double **grid, int *grid_size, int num, int k) ;

void grid_extend(double x, double **grid, int *grid_size, int *num, int k, double **coef) ;


// 反向传播更新 coef
void backward_neural(Neural *neural, double x, double error, double learning_rate, double lambda) ;

void calculate_previous_error(Neural* neural, double current_error) ;

void fix_fun(Neural *neural, double *inputs, double *outputs, int num_samples) ;
// 封装的训练函数
void train(Neural *neural, double *inputs, double *outputs, int num_samples, int num_epochs, double learning_rate, double lambda) ;

void compute_loss(Neural *neural, double *inputs, double *outputs, int num_samples) ;

#endif