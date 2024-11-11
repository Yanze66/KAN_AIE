#include <stdio.h>
#include <stdlib.h>
#include <math.h>



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

//  define the function to calculate the B-spline basis function,x为输入，grid为节点，num为节点数，k为阶数，result为输出
void B_batch(double x, double *grid, int num, int k, double *result) {
    int num_basis = num + k - 1;
    double x_val = x;
    for (int i = 0; i < num_basis; i++) {
        result[i] = calculate_bspline(grid, k, i, x_val);
    }
}


// 初始化神经元以便使用样条
typedef struct {
    double *coef;
    double *grid;
    int num; // 样条节点数
    int k;   // 样条阶数

} Neural;

Neural* init_neural(int num, int k) {
    Neural* neural = (Neural*)malloc(sizeof(Neural));
    neural->coef = (double*)malloc((num + k - 1) * sizeof(double));
    neural->grid = (double*)malloc((num + 2 * k) * sizeof(double));
    neural->num = num;
    neural->k = k;

    // 初始化系数为1（测试用）
    for (int i = 0; i < num + k - 1; i++) {
        neural->coef[i] = 0.0;
    }

    // 初始化网格
    double step = 2.0 / (num + 1); // 在[-1, 1]之间均匀分布
    for (int i = 0; i < num + 2 * k; i++) {
        neural->grid[i] = -1.0 + i * step;
        printf("%f\n", neural->grid[i]);
    }

    return neural;
}

double forward_neural(Neural *neural, double x) {
    double sx_result[neural->num + neural->k];
    B_batch(x, neural->grid, neural->num, neural->k, sx_result);
    double result = 0.0;
    for (int i = 0; i < neural->num+neural->k-1; i++) {
        result += neural->coef[i] * sx_result[i];
    }
    return result;
}

// 反向传播更新 coef
void backward_neural(Neural *neural, double x, double error, double learning_rate) {
    double sx_result[neural->num + neural->k];
    B_batch(x, neural->grid, neural->num, neural->k, sx_result);

    for (int i = 0; i < neural->num +neural->k-1; i++) {
        double gradient = error * sx_result[i];
        neural->coef[i] -= learning_rate * gradient;
    }
}

int main() {
    int num = 150;
    int k = 3;
    int num_samples = 100;
    int num_epochs = 1000;
    double learning_rate = 0.01;

    Neural *neural = init_neural(num, k);

    for (int epoch = 0; epoch < num_epochs; epoch++) {
        double total_loss = 0.0;
        double adjusted_learning_rate = learning_rate / (1 + 0.001 * epoch); // 衰减学习率
        for (int i = 0; i < num_samples; i++) {
            double x = ((double)rand() / RAND_MAX) * 2 - 1; // 随机 x
            double true_value = x*x; // 真实函数值

            double predicted_value = forward_neural(neural, x);
            double error = predicted_value - true_value;

            total_loss += error * error;

            // 执行反向传播，更新 coef
            backward_neural(neural, x, error, adjusted_learning_rate);
        }

        total_loss = sqrt(total_loss / num_samples);
        printf("Epoch %d: RMSE Loss = %f\n", epoch + 1, total_loss);
    }

    free(neural->coef);
    free(neural->grid);
    free(neural);

    return 0;
}