#include <stdio.h>
#include <time.h>
#include <x86intrin.h>
// #include "../spline2fun_adam.h"
#include "../spline2fun_loop.h"
// 声明 B_batch 函数
void test_b_batch() {
    // Define grid dimensions
    int num = 5;
    int k = 3;
    int num_interval = num + 2 * k;
    int num_basis = num_interval - k - 1;
    double test_values[] = {-0.5, 0.0, 0.5};
    int num_tests = sizeof(test_values) / sizeof(test_values[0]);

    // Allocate 2D grid array (for each test value) and result array
    double **grid_array = (double **)malloc(num_tests * sizeof(double *));
    double **result_array = (double **)malloc(num_tests * sizeof(double *));
    for (int j = 0; j < num_tests; j++) {
        grid_array[j] = (double *)malloc(num_interval * sizeof(double));
        result_array[j] = (double *)malloc(num_basis * sizeof(double));

        // Fill grid for each test case with uniformly distributed points in [-1, 1]
        for (int i = 0; i < num_interval; i++) {
            grid_array[j][i] = -1.0 + i * (2.0 / (num_interval - 1));
            printf("grid_array[%d][%d] = %f\n", j, i, grid_array[j][i]);
        }
    }

    // Call B_batch with test values
    B_batch(test_values, num_tests, grid_array, num, 4, result_array);

    // Print results for each test value
    for (int j = 0; j < num_tests; j++) {
        printf("B-spline basis functions for x = %.2f:\n", test_values[j]);
        for (int i = 0; i < num_basis; i++) {
            printf("B[%d] = %f\n", i, result_array[j][i]);
        }
        printf("\n");
    }

    // Free allocated memory
    // for (int j = 0; j < num_tests; j++) {
    //     free(grid_array[j]);
    //     free(result_array[j]);
    // }
    free(grid_array);
    free(result_array);
}

void test_grid_extend() {
    int k = 3; // 样条阶数
    int initial_num = 5; // 初始节点数量
    int grid_size =  initial_num + 2*k; // 初始网格大小
    double *grid = (double *)malloc(grid_size * sizeof(double));
    double *coef = (double *)calloc(initial_num + k - 1, sizeof(double));

    // 初始化网格为 [-1.0, 1.0] 区间内均匀分布的点
    for (int i = 0; i < grid_size; i++) {
        grid[i] = -1.0 + i * (2.0 / (grid_size - 1));
    }

    // 打印初始网格
    printf("Initial grid:\n");
    for (int i = 0; i < grid_size; i++) {
        printf("%f ", grid[i]);
    }
    printf("\n");

    // 测试用例1: x值在网格左侧超出
    double x1 = -1.5;
    grid_extend(x1, &grid, &grid_size, &initial_num, k, &coef);
    printf("Extended grid for x = %.2f:\n", x1);
    for (int i = 0; i < grid_size; i++) {
        printf("%f ", grid[i]);
    }
    printf("\n");

    // 测试用例2: x值在网格右侧超出
    double x2 = 1.5;
    grid_extend(x2, &grid, &grid_size, &initial_num, k, &coef);
    printf("Extended grid for x = %.2f:\n", x2);
    for (int i = 0; i < grid_size; i++) {
        printf("%f ", grid[i]);
    }
    printf("\n");

    // 测试用例3: x值在网格内
    double x3 = 0.5;
    grid_extend(x3, &grid, &grid_size, &initial_num, k, &coef);
    printf("Extended grid for x = %.2f:\n", x3);
    for (int i = 0; i < grid_size; i++) {
        printf("%f ", grid[i]);
    }
    printf("\n");

    printf("initial_num = %d\n", initial_num);

    // 释放内存
    free(grid);
    free(coef);
}

int main() {
    // 定义变量用于存储时间和CPU周期
    clock_t start_time, end_time;
    unsigned long long start_cycles, end_cycles;

    // 记录开始时间和CPU周期
    start_time = clock();
    start_cycles = __rdtsc();

    // 运行测试函数
    test_b_batch();
    // test_grid_extend();

    // 记录结束时间和CPU周期
    end_cycles = __rdtsc();
    end_time = clock();

    // 计算并输出时间和周期
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", time_spent);
    printf("CPU cycles: %llu\n", end_cycles - start_cycles);

    return 0;
}