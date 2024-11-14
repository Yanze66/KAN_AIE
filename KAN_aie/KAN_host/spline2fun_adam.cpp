#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "spline2fun_adam.h"
// #include "kanlayer.h"

// double calculate_bspline(double *grid, int num_of_basis, double x) {
//     double N0[5] = {0};
//     double N1[4] = {0};
//     double N2[3] = {0};
//     double N3[2] = {0};
//     double N4[1] = {0};

//     // 计算0阶B样条系数
//     for (int i = 0; i <= 4; i++) {
//         N0[i] = (grid[i] <= x && x < grid[i + 1]) ? 1.0 : 0.0;
//     }

//     // 计算1阶B样条系数
//     for (int i = 0; i < 4; i++) {
//         double alpha_denominator = grid[i + 1] - grid[i];
//         double beta_denominator = grid[i + 2] - grid[i + 1];
//         double alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N0[i];
//         double beta = (beta_denominator == 0) ? 0 : (grid[i + 2] - x) / beta_denominator * N0[i + 1];
//         N1[i] = alpha + beta;
//     }

//     // 计算2阶B样条系数
//     for (int i = 0; i < 3; i++) {
//         double alpha_denominator = grid[i + 2] - grid[i];
//         double beta_denominator = grid[i + 3] - grid[i + 1];
//         double alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N1[i];
//         double beta = (beta_denominator == 0) ? 0 : (grid[i + 3] - x) / beta_denominator * N1[i + 1];
//         N2[i] = alpha + beta;
//     }

//     // 计算3阶B样条系数
//     for (int i = 0; i < 2; i++) {
//         double alpha_denominator = grid[i + 3] - grid[i];
//         double beta_denominator = grid[i + 4] - grid[i + 1];
//         double alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N2[i];
//         double beta = (beta_denominator == 0) ? 0 : (grid[i + 4] - x) / beta_denominator * N2[i + 1];
//         N3[i] = alpha + beta;
//     }

//     // 计算4阶B样条系数
//     for (int i = 0; i < 1; i++) {
//         double alpha_denominator = grid[i + 4] - grid[i];
//         double beta_denominator = grid[i + 5] - grid[i + 1];
//         double alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N3[i];
//         double beta = (beta_denominator == 0) ? 0 : (grid[i + 5] - x) / beta_denominator * N3[i + 1];
//         N4[i] = alpha + beta;
//     }

//     return N4[0];
// }

// void B_batch(double x, double *grid, int num, int k, double *result) {
//     int num_basis = num + k - 1;
//     double x_val = x;
//     for (int i = 0; i < num_basis; i++) {
//         result[i] = calculate_bspline(grid + i, 4, x_val);
//     }
// }


//  define the function to calculate the B-spline basis function
double calculate_bspline(double *grid, int k, int num_of_basis, double x) {
    if(k==0){
        return (grid[num_of_basis] <= x && x < grid[num_of_basis + 1]) ? 1.0 : 0.0;
    } else {
        double alpha_denominator = grid[num_of_basis + k] - grid[num_of_basis];
        double beta_denominator = grid[num_of_basis + k + 1] - grid[num_of_basis + 1];
        double alpha = (alpha_denominator == 0) ? 0 : (x - grid[num_of_basis]) / alpha_denominator * calculate_bspline(grid, k - 1, num_of_basis, x);
        double beta = (beta_denominator == 0) ? 0 : (grid[num_of_basis + k + 1] - x) / beta_denominator * calculate_bspline(grid, k - 1, num_of_basis + 1, x);
        return alpha + beta;
    }
}

//  define the function to calculate the B-spline basis function,x为输入，grid为节点，num_interval为节点数，k为阶数，result为输出
void B_batch(double x, double *grid, int num_interval, int k, double *result) {
    int num_basis = num_interval + k - 1;
    double x_val = x;
    for (int i = 0; i < num_basis; i++) {
        result[i] = calculate_bspline(grid, k, i, x_val);
    }
}






double calculate_bspline_derivative(double *grid, int k, int num_of_basis, double x) {
    if (k == 0) {
        return 0.0; // 因为0阶B样条是常数，导数为0
    } else {
        double alpha_denominator = grid[num_of_basis + k] - grid[num_of_basis];
        double beta_denominator = grid[num_of_basis + k + 1] - grid[num_of_basis + 1];

        double alpha_term = 0, beta_term = 0;
        if (alpha_denominator != 0) {
            alpha_term = (k / alpha_denominator) * calculate_bspline(grid, k - 1, num_of_basis, x);
        }
        
        if (beta_denominator != 0) {
            beta_term = (-k / beta_denominator) * calculate_bspline(grid, k - 1, num_of_basis + 1, x);
        }

        return alpha_term + beta_term;
    }
}

double relu(double x) {
    return x > 0 ? x : 0;
}

double silu(double x) {
    return x / (1 + exp(-x));
}

double silu_derivative(double x) {
    double sigmoid = 1 / (1 + exp(-x));
    return sigmoid + x * sigmoid * (1 - sigmoid);
}

// Extend the grid to include `x` if it's outside the current range

void grid_extend(double x, double **grid, int *grid_size, int *num, int k, double **coef) {
    double min_grid = (*grid)[0];
    double max_grid = (*grid)[*grid_size - 1];
    double current_step = (max_grid - min_grid) / (*grid_size - 1);
    // 当 x 超出当前网格边界时
    if (x < min_grid || x > max_grid) {

    double new_min = min_grid;
    double new_max = max_grid;

        // 如果 x 小于 min_grid，需要向左扩展
        if (x < min_grid) {
            new_min = x - current_step;  // 确保最小值覆盖到 x
        }

        // 如果 x 大于 max_grid，需要向右扩展
        if (x > max_grid) {
            new_max = x + current_step;  // 确保最大值超过 x
        }

        // 计算新的节点数量
        int additional_nodes_left = (int)ceil((min_grid - new_min) / current_step);
        int additional_nodes_right = (int)ceil((new_max - max_grid) / current_step);
        int new_num = *num + additional_nodes_left + additional_nodes_right;

        // 重新计算新的网格大小
        int new_grid_size = new_num + 2* k;
        int new_coef_size = new_num + k - 1;

        // 为新网格和系数分配内存
        double *new_grid = (double *)malloc(new_grid_size * sizeof(double));
        double *new_coef = (double *)calloc(new_coef_size, sizeof(double));

        // 填充新网格
        for (int i = 0; i < new_grid_size; i++) {
            new_grid[i] = new_min + i * current_step;
        }
       

        
        for (int i = 0; i < new_coef_size; i++) {
            new_coef[i] = normal_random(0.1);
        }
        

        // 释放旧的内存
        free(*grid);
        free(*coef);
        
        *grid = new_grid;
        *coef = new_coef;
        *grid_size = new_grid_size;
        *num = new_num; // 更新 num
    }
}



// define the structure of a single neuron
Neural* init_neural(int l,int in_index, int out_index, int num, int k) {
    Neural* neural = (Neural*)malloc(sizeof(Neural));
    neural->layer = l;
    neural->in_index = in_index;
    neural->out_index = out_index;
  
    neural->history_index = 0;
    neural->output_history = (double*)calloc(20, sizeof(double));
    for(int i = 0; i < 20; i++){
        neural->output_history[i] = 0.0;
    }

    neural->coef = (double*)malloc((num + k - 1) * sizeof(double));
    neural->grid = (double*)malloc((num + 2 * k) * sizeof(double));
    neural->num = num;
    neural->k = k;
    neural->grid_size = num + 2 * k;
    
    
    neural->wb = ((double)rand() / RAND_MAX) ;

    // initialize the coefficients to 1 (for testing)
    for (int i = 0; i < num + k - 1; i++) {
        neural->coef[i] = normal_random(0.1);
    }

    // initialize the grid
    double step = 10.0 / (num + 2*k); // evenly distributed in [-1, 1]
    for (int i = 0; i < num + 2 * k; i++) {
        //neural->grid[i] = -10.0 + i * step;
        neural->grid[i] = -5 + i * step;
    }

    neural->in_value = 0.0;
    neural->out_value = 0.0;
    neural->active = 1;

    neural->fun_fixed = 0;
    neural->fx = NULL;
    neural->fx_derivative = NULL;
    neural->best_a = 0;
    neural->best_b = 0;
    neural->best_c = 0;
    neural->best_d = 0; 

    neural->pre_error = 0.0;
    return neural;
}

double forward_neural(Neural *neural, double x) {
    double result = 0.0;
    neural->in_value = x;
    grid_extend(x, &neural->grid, &neural->grid_size, &neural->num, neural->k, &neural->coef);
    if(neural->active == 0){
        return 0.0;
    }else {
        if(!neural->fun_fixed){
            double sx_result[neural->num + neural->k];
            B_batch(x, neural->grid, neural->num, neural->k, sx_result);
            double spline_value = 0.0;
    
            for (int i = 0; i < neural->num+neural->k-1; i++) {
                spline_value += neural->coef[i] * sx_result[i];
            }
            //result = neural->active * (relu(x) * neural->wb + spline_value);
            result = neural->active * (silu(x) * neural->wb + spline_value);
        }else if(neural->fun_fixed){
            result = neural->fx(neural->best_a * x + neural->best_b) * neural->best_c + neural->best_d;
        }
        neural->out_value = result;
        return result;
    }
}


// BP update coef, SG version, L1 regularization
void backward_neural(Neural *neural, double x, double error, double learning_rate, double lambda) {
    if(neural->active == 0){
        return;
    }else{
        if(!neural->fun_fixed){
        double sx_result[neural->num + neural->k-1];
        B_batch(x, neural->grid, neural->num, neural->k, sx_result);
        
        
        for (int i = 0; i < neural->num +neural->k-1; i++) {
            double gradient = error * sx_result[i];
            // L1正则化梯度调整：在原有梯度上加上 L1 的部分
            if (neural->coef[i] > 0) {
                gradient += lambda;
            } else if (neural->coef[i] < 0) {
                gradient -= lambda;
            }

            neural->coef[i] -= learning_rate * gradient;
        }

         // 更新权重 wb
        //double activation_derivative = (x > 0) ? 1 : 0; // ReLU derivative
        double activation_derivative = silu_derivative(x); // silu derivative
        double gradient_wb = error * activation_derivative * x ;
         // 对偏置项使用 L1 正则化
        if (neural->wb > 0) {
            gradient_wb += lambda;
        } else if (neural->wb < 0) {
            gradient_wb -= lambda;
        }

        neural->wb -= learning_rate * gradient_wb;
        return;  
        } else if(neural->fun_fixed){
            
        double fx_derivative = neural->fx_derivative(neural->best_a * x + neural->best_b);
        double gradient_a = error * neural->best_c* fx_derivative * x;
        double gradient_b = error * neural->best_c*fx_derivative;
        double gradient_c = error * fx_derivative;
        double gradient_d = error;

        // // 添加 L2 正则化项
        // gradient_a += lambda * neural->best_a;
        // gradient_b += lambda * neural->best_b;
        // gradient_c += lambda * neural->best_c;
        // gradient_d += lambda * neural->best_d;

        // 应用梯度裁剪
        double max_grad = 1.0;
        if (gradient_a > max_grad) gradient_a = max_grad;
        if (gradient_a < -max_grad) gradient_a = -max_grad;
        if (gradient_b > max_grad) gradient_b = max_grad;
        if (gradient_b < -max_grad) gradient_b = -max_grad;
        if (gradient_c > max_grad) gradient_c = max_grad;
        if (gradient_c < -max_grad) gradient_c = -max_grad;
        if (gradient_d > max_grad) gradient_d = max_grad;
        if (gradient_d < -max_grad) gradient_d = -max_grad;

        double new_learning_rate = learning_rate /500;
        // 应用梯度更新
        neural->best_a -= new_learning_rate * gradient_a;
        neural->best_b -= new_learning_rate * gradient_b;
        neural->best_c -= new_learning_rate * gradient_c;
        neural->best_d -= new_learning_rate * gradient_d;

    }
    
    }
    
}

// Calculate the previous error for the neuron

void calculate_previous_error(Neural* neural, double current_error) {
    if (neural->active == 0) {
        neural->pre_error = 0;
        return;
    } else if(!neural->fun_fixed){
        double sigmoid = 1.0 / (1.0 + exp(-neural->in_value));
        double silu_derivative = sigmoid + (neural->in_value * sigmoid * (1.0 - sigmoid));
        
        // Calculate B-Spline derivatives at x = neural->in_value
        double sx_derivative[neural->num + neural->k];
        for (int i = 0; i < neural->num + neural->k - 1; i++) {
            sx_derivative[i] = calculate_bspline_derivative(neural->grid, neural->k, i, neural->in_value);
        }

        // Combine derivatives to calculate the previous error
        double spline_derivative = 0.0;
        for (int i = 0; i < neural->num + neural->k - 1; i++) {
            spline_derivative += neural->coef[i] * sx_derivative[i];
        }

        double previous_error = current_error * (neural->wb * silu_derivative + spline_derivative);
        neural->pre_error = previous_error;
        return;
    }else{
        double fx_derivative = neural->fx_derivative(neural->best_a * neural->in_value + neural->best_b);
        double previous_error = current_error * neural->best_c * fx_derivative;
        neural->pre_error = previous_error;
        return;
    }
}

// 修正函数 y = c f(ax + b) + d
  
    // 定义可能的函数选择
    //function[0] = identity
    //function[1] = square
    //function[2] = cube
    //function[3] = power4
    //function[4] = power5
    //function[5] = reciprocal
    //function[6] = reciprocal_square
    //function[7] = reciprocal_cube
    //function[8] = reciprocal_x4
    //function[9] = reciprocal_x5
    //function[10] = reciprocal_sqrt
    //function[11] = exp_func
    //function[12] = log_func
    //function[13] = abs_func
    //function[14] = sin_func
    //function[15] = cos_func
    //function[16] = tanh_func
    //function[17] = sgn_func
    //function[18] = arctan_func
    //function[19] = sqrt_func
    //function[20] = cube_root
    //function[21] = arctanh_func
    //function[22] = arccos_func
    //function[23] = arcsin_func
    //function[24] = zero_func
    //function[25] = gauss_func
    //function[26] = tan_func
void fix_fun(Neural *neural, double *inputs, double *outputs, int num_samples) {
    if(neural->active == 0){
        return;
    }
    
    //double best_a = 0, best_b = 0, best_c = 0, best_d = 0;
    double min_rmse = INFINITY;
    double a_range = 6, b_range = 6, step_ab = 0.1,d_range = 2,step_cd = 0.5;
  

    double (*functions[])(double) = {  identity, square, cube, power4, power5, reciprocal, reciprocal_square, reciprocal_cube, reciprocal_x4, reciprocal_x5,  reciprocal_sqrt, 
    exp_func, log_func, abs_func, sin_func, cos_func, tanh_func, sgn_func, arctan_func, sqrt_func, cube_root, arctanh_func,arccos_func, arcsin_func, zero_func, gauss_func,tan_func};
    double (*functions_derivative[])(double) = {  identity_derivative, square_derivative, cube_derivative, power4_derivative, power5_derivative, reciprocal_derivative, reciprocal_square_derivative,
     reciprocal_cube_derivative, reciprocal_x4_derivative, reciprocal_x5_derivative,  reciprocal_sqrt_derivative, exp_func_derivative, log_func_derivative, abs_func_derivative, sin_func_derivative, cos_func_derivative, 
     tanh_func_derivative, sgn_func_derivative, arctan_func_derivative, sqrt_func_derivative, cube_root_derivative, arctanh_func_derivative,arccos_func_derivative, arcsin_func_derivative, zero_func_derivative, gauss_func_derivative}; 
    int num_functions = sizeof(functions) / sizeof(functions[0]);
    int target_function = 0;
    for (int f = 0; f < num_functions; f++) {
        for (double a = -a_range; a <= a_range; a += step_ab) {
            for (double b = -a_range; b <= a_range; b += step_ab) {
                for (double c = -d_range; c <= d_range; c += step_cd) {
                    for (double d = -d_range; d <= d_range; d += step_cd) {
                        double total_error = 0.0;
                        for (int i = 0; i < num_samples; i++) {
                            double scaled_input = a * inputs[i] + b;
                            double output_pred = c * functions[f](scaled_input) + d;
                            //double output_pred = functions[f](scaled_input);
                            double error = outputs[i] - output_pred;
                            total_error += error * error;
                        }
                        double rmse = sqrt(total_error / num_samples);
                        if (rmse < min_rmse) {
                            min_rmse = rmse;
                            neural->best_a = a;
                            neural->best_b = b;
                            neural->best_c = c;
                            neural->best_d = d;
                            neural->fx = functions[f];
                            neural->fx_derivative = functions_derivative[f];
                            target_function = f;
                            //printf("a = %.2f, b = %.2f, c = %.2f, d = %.2f, function[%d] = %p, rmse = %.2f\n", a, b, c, d,f, functions[f], rmse);
                        }
                    }
                }
            }
        }
    }

    printf("Optimal parameters found: a = %.2f, b = %.2f, c = %.2f, d = %.2f,function[%d]= %p,rmse = %.2f\n", neural->best_a, neural->best_b, neural->best_c, neural->best_d,target_function,neural->fx, min_rmse);

    // set fun_fixed to 1 to indicate that the function is fixed
    neural->fun_fixed = 1;
}

// train the neural network
void train(Neural *neural, double *inputs, double *outputs, int num_samples, int num_epochs, double learning_rate, double lambda) {
    for (int epoch = 0; epoch < num_epochs; epoch++) {
        double total_loss = 0.0;
        double adjusted_learning_rate = learning_rate / (1 + 0.001 * epoch);
        //double adjusted_learning_rate = learning_rate;
        for (int i = 0; i < num_samples; i++) {
            double x = inputs[i];
            double true_value = outputs[i];

            double predicted_value = forward_neural(neural, x);
            
            //printf("Sample %d, x = %.3f, real-output=%.3f, predicted=%.3f\n", i, x, true_value, predicted_value);
            

            double error = predicted_value - true_value;
            
            total_loss += error * error;

            backward_neural(neural, x, error, adjusted_learning_rate,lambda);

        }
        
        total_loss = sqrt(total_loss / num_samples);
        printf("Epoch %d: RMSE Loss = %f\n", epoch + 1, total_loss);
    }
}

void compute_loss(Neural *neural, double *inputs, double *outputs, int num_samples) {
    double total_loss = 0;

    for (int i = 0; i < num_samples; i++) {
        double predicted_value = forward_neural(neural, inputs[i]);
        double error = predicted_value - outputs[i];
        total_loss += error * error;
    }

    total_loss = sqrt(total_loss / num_samples);
    printf("Final RMSE Loss: %f\n", total_loss);
}

