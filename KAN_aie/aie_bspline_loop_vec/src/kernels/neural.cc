// #include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <string.h>
#include <cstdio>
#include <stdio.h>
// #include <math.h>
// #include <time.h>
#include <adf.h>
#include <aie_api/aie.hpp>
#include "aie_api/aie_types.hpp"
#include "../kernels.h"
// #include "adf/x86sim/streamApi.h"
// #include "adf/x86sim/streamApi.h"

#define ctl_num  5
#define order 3




float calculate_bspline(float *grid, int num_of_basis, float x) {
    float N0[5] = {0};
    float N1[4] = {0};
    float N2[3] = {0};
    float N3[2] = {0};
    float N4[1] = {0};

    // 计算0阶B样条系数
    for (int i = 0; i <= 4; i++) {
        N0[i] = (grid[i] <= x && x < grid[i + 1]) ? 1.0 : 0.0;
    }

    // 计算1阶B样条系数
    for (int i = 0; i < 4; i++) {
        float alpha_denominator = grid[i + 1] - grid[i];
        float beta_denominator = grid[i + 2] - grid[i + 1];
        float alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N0[i];
        float beta = (beta_denominator == 0) ? 0 : (grid[i + 2] - x) / beta_denominator * N0[i + 1];
        N1[i] = alpha + beta;
    }

    // 计算2阶B样条系数
    for (int i = 0; i < 3; i++) {
        float alpha_denominator = grid[i + 2] - grid[i];
        float beta_denominator = grid[i + 3] - grid[i + 1];
        float alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N1[i];
        float beta = (beta_denominator == 0) ? 0 : (grid[i + 3] - x) / beta_denominator * N1[i + 1];
        N2[i] = alpha + beta;
    }

    // 计算3阶B样条系数
    for (int i = 0; i < 2; i++) {
        float alpha_denominator = grid[i + 3] - grid[i];
        float beta_denominator = grid[i + 4] - grid[i + 1];
        float alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N2[i];
        float beta = (beta_denominator == 0) ? 0 : (grid[i + 4] - x) / beta_denominator * N2[i + 1];
        N3[i] = alpha + beta;
    }

    // 计算4阶B样条系数
    for (int i = 0; i < 1; i++) {
        float alpha_denominator = grid[i + 4] - grid[i];
        float beta_denominator = grid[i + 5] - grid[i + 1];
        float alpha = (alpha_denominator == 0) ? 0 : (x - grid[i]) / alpha_denominator * N3[i];
        float beta = (beta_denominator == 0) ? 0 : (grid[i + 5] - x) / beta_denominator * N3[i + 1];
        N4[i] = alpha + beta;
    }

    return N4[0];
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
            alpha_term = (k / alpha_denominator) * calculate_bspline(grid,  num_of_basis, x);
        }   
        if (beta_denominator != 0) {
            beta_term = (-k / beta_denominator) * calculate_bspline(grid,  num_of_basis + 1, x);
        }
        return alpha_term + beta_term;
    }
}

inline void B_batch(float x, float * __restrict grid, uint32_t num, uint32_t k, float * __restrict result ){
    uint32_t num_basis = num + k - 1;
    float x_val = x;
    for (int i = 0; i < num_basis; i++) {
        result[i] = calculate_bspline(grid+i,  4, x_val);
    }    
}


/*************************************************/

// GMIO is suppoused to transfer the data in multiple of 32bit, or it renders deadlock when reading.
void spline1(input_stream<float> *  __restrict datain, output_stream<float>*  __restrict dataout,output_stream<float>*  __restrict resultout)
{
    
    char pipe =4;
    char num =5;
    char k = 3;
    
    int num_interval = num+ 2*k;
    float grid[num_interval];
    float result[num_interval - k - 1]; // num_basis = num_interval - k - 1

    // 填充网格，在 [-1, 1] 间均匀填充 11 个点
    for (int i = 0; i < num_interval; i++) {
        grid[i] = -1.0 + i * (2.0 / (num_interval - 1));
    }

    //read data
    float buf[pipe];
    for(int i=0;i<pipe;i++){
        buf[i]=readincr(datain);
    }
    
    //pickup input
    float x = buf[0];

    //transport vector
    for(int i=0;i<pipe;i++){
        writeincr(dataout, buf[i]);
    }

    //calculate
    B_batch(x,grid,num,k,result);//order = k
    


    //transport result
    
        // printf("\ninput %d\t",i);
        // for(int j=0;j<num+k-1;j++){           
        //     printf("value %f,",result[j]);
        //     writeincr(bufout, result[j]);
        // }
    
    
    //output result result
    for(int j=0;j<num+k-1;j++){           
            printf("value %f,",result[j]);
            writeincr(resultout, result[j]);
        }

    writeincr(resultout,1010); //fill the buf in length of 128-bit

}

void spline2(input_stream<float> *  __restrict datain, output_stream<float>*  __restrict dataout,input_stream<float>*  __restrict resultin, output_stream<float>*  __restrict resultout)
{
    
    char pipe =4;
    char num =5;
    char k = 3;
    
    int num_interval = num+ 2*k;
    float grid[num_interval];
    float result[num_interval - k - 1]; // num_basis = num_interval - k - 1

    // 填充网格，在 [-1, 1] 间均匀填充 11 个点
    for (int i = 0; i < num_interval; i++) {
        grid[i] = -1.0 + i * (2.0 / (num_interval - 1));
    }

    //read data
    float buf[pipe];
    for(int i=0;i<pipe;i++){
        buf[i]=readincr(datain);
    }
    
    //pickup input
    float x = buf[1];

    //transport vector
    for(int i=0;i<pipe;i++){
        writeincr(dataout, buf[i]);
    }

    //calculate
    B_batch(x,grid,num,k,result);//order = k
    


    //transport result
    float result_reg;
    for(int j=0;j<1;j++){  //only accept 1 level
        for(int j=0;j<num+k;j++){           
            // printf("value %f,",result[j]);
            result_reg = readincr(resultin);
            writeincr(resultout, result_reg);
        }
    }
        
    
    
    //output result result
    for(int j=0;j<num+k-1;j++){           
            printf("value %f,",result[j]);
            writeincr(resultout, result[j]);
        }

    writeincr(resultout,1010); //fill the buf in length of 128-bit

}

void spline3(input_stream<float> *  __restrict datain, input_stream<float>*  __restrict resultin, output_stream<float>*  __restrict resultout)
{
    
    char pipe =4;
    char num =5;
    char k = 3;
    
    int num_interval = num+ 2*k;
    float grid[num_interval];
    float result[num_interval - k - 1]; // num_basis = num_interval - k - 1

    // 填充网格，在 [-1, 1] 间均匀填充 11 个点
    for (int i = 0; i < num_interval; i++) {
        grid[i] = -1.0 + i * (2.0 / (num_interval - 1));
    }

    //read data
    float buf[pipe];
    for(int i=0;i<pipe;i++){
        buf[i]=readincr(datain);
    }
    
    //pickup input
    float x = buf[2];

    //transport vector
    // for(int i=0;i<pipe;i++){
    //     writeincr(dataout, buf[i]);
    // }

    //calculate
    B_batch(x,grid,num,k,result);//order = k
    


    //transport result
    float result_reg;
    for(int j=0;j<2;j++){  // accept 2 level
        for(int j=0;j<num+k;j++){           
            // printf("value %f,",result[j]);
            result_reg = readincr(resultin);
            writeincr(resultout, result_reg);
        }
    }
        
    
    
    //output result result
    for(int j=0;j<num+k-1;j++){           
            printf("value %f,",result[j]);
            writeincr(resultout, result[j]);
        }

    writeincr(resultout,1010); //fill the buf in length of 128-bit

}