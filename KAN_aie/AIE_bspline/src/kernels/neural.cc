// #include <stdio.h>
// #include <stdlib.h>
// #include <cstdint>
// #include <string.h>
// #include <cstdio>
// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <adf.h>
// #include <aie_api/aie.hpp>
// #include "aie_api/aie_types.hpp"
#include "../kernels.h"
// #include "adf/x86sim/streamApi.h"

#define ctl_num  5
#define order 3




inline float calculate_bspline(float *__restrict grid, uint32_t k, uint32_t num_of_basis, float x) {
    if(k==0){
        return (grid[num_of_basis] <= x && x < grid[num_of_basis + 1]) ? 1.0 : 0.0;
    } else {
        float alpha_denominator = grid[num_of_basis + k] - grid[num_of_basis];
        float beta_denominator = grid[num_of_basis + k + 1] - grid[num_of_basis + 1];
        float alpha = (alpha_denominator == 0) ? 0 : (x - grid[num_of_basis]) / alpha_denominator * calculate_bspline(grid, k - 1, num_of_basis, x);
        float beta = (beta_denominator == 0) ? 0 : (grid[num_of_basis + k + 1] - x) / beta_denominator * calculate_bspline(grid, k - 1, num_of_basis + 1, x);
        return alpha + beta;
    }
}

inline float calculate_bspline_derivative(float *__restrict grid, uint32_t k, uint32_t num_of_basis, float x) {
    if (k == 0) {
        return 0.0; // derivative is 0 for constant
    } else {
        float alpha_denominator = grid[num_of_basis + k] - grid[num_of_basis];
        float beta_denominator = grid[num_of_basis + k + 1] - grid[num_of_basis + 1];

        float alpha_term = 0;
        float beta_term = 0;
        if (alpha_denominator != 0) {
            alpha_term = (k / alpha_denominator) * calculate_bspline(grid, k - 1, num_of_basis, x);
        }   
        if (beta_denominator != 0) {
            beta_term = (-k / beta_denominator) * calculate_bspline(grid, k - 1, num_of_basis + 1, x);
        }
        return alpha_term + beta_term;
    }
}

inline void B_batch(float x, float * __restrict grid, uint32_t num, uint32_t k, float * __restrict result ){
    uint32_t num_basis = num + k - 1;
    float x_val = x;
    for (int i = 0; i < num_basis; i++) {
        result[i] = calculate_bspline(grid, k, i, x_val);
    }    
}


/*************************************************/

// GMIO is suppoused to transfer the data in multiple of 32bit, or it renders deadlock when reading.
void spline(input_stream<float> *  __restrict bufin, output_stream<float>*  __restrict bufout)
{
    

    char num =5;
    char k = 3;
    
    int num_interval = num+ 2*k;
    float grid[num_interval];
    float result[num_interval - k - 1]; // num_basis = num_interval - k - 1

    // 填充网格，在 [-1, 1] 间均匀填充 11 个点
    for (int i = 0; i < num_interval; i++) {
        grid[i] = -1.0 + i * (2.0 / (num_interval - 1));
    }

    
    printf("start aie");

    
    float z = readincr(bufin);
    printf("\nz=%f\n",z);
   

    //float result[num+k-1];
    for(int i=0;i<3;i++){
        printf("start loop");
        float x = readincr(bufin);
        printf("\nx=%f\n",x);        
        B_batch(x,grid,num,4,result);
        // printf("\ninput %d\t",i);
        for(int j=0;j<num+k-1;j++){           
            printf("value %f,",result[j]);
            writeincr(bufout, result[j]);
        }
    }
    
    



}