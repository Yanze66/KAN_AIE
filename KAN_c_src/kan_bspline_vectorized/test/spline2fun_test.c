#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../utils.h"
#include "../spline2fun_adam.h"
// #include "../spline2fun.h"

int main() {
    // Initialize random seed
    srand(time(NULL));

    int num = 50;
    int k = 3;
    int num_samples = 3000;
    int num_epochs = 100;
    double learning_rate = 0.01;
    double lambda = 0.001; // regularization parameter


    Neural *neural = init_neural(0,1,1,num, k);

    // 设置生成的数据集
    double inputs[num_samples];
    double outputs[num_samples];
    generate_data(inputs, outputs, num_samples);

    //standardize_data(inputs, num_samples);
    //standardize_data(outputs, num_samples);
    

    // 执行训练过程
    train(neural, inputs, outputs, num_samples, num_epochs, learning_rate,lambda);
    


    // 测试训练后的网络
    double sample_inputs[10];
    double sample_outputs[10];
    double predicted_out[10];
    generate_data(sample_inputs, sample_outputs, 10);

    for(int i=0;i<10;i++){       
        predicted_out[i] = forward_neural(neural, sample_inputs[i]);
        printf("Sample %d, x = %.3f, result =%.3f, Predicted = %.3f\n", i, sample_inputs[i], sample_outputs[i],predicted_out[i]);
    }
    compute_loss(neural, sample_inputs, sample_outputs, 10);




    //fixing neural network 
    printf("begins to fix function: \n");
    double fix_inputs[20];
    double fix_outputs[20];
    generate_data(fix_inputs, fix_outputs, 20);
    
    for(int i=0;i<20;i++){
        fix_outputs[i] = forward_neural(neural, fix_inputs[i]);
    }
    fix_fun(neural, fix_inputs, fix_outputs, 20);

    printf("Finishing fixing, begins to compute loss: \n");

    //test the fixed function
    generate_data(sample_inputs, sample_outputs, 10);
    //standardize_data(sample_inputs, 20);
    //standardize_data(sample_outputs, 20);
    for(int i=0;i<10;i++){       
        predicted_out[i] = forward_neural(neural, sample_inputs[i]);
        printf("Sample %d, x = %.3f, result =%.3f, Predicted = %.3f\n", i, sample_inputs[i], sample_outputs[i],predicted_out[i]);
    }
    
    compute_loss(neural, sample_inputs, sample_outputs, 10);


    //training the fixed function
    
    // double inputs_after_fix[1000];
    // double outputs_after_fix[1000];
    // generate_data(inputs_after_fix, outputs_after_fix, 1000);
    // // 执行训练过程
    // train(neural, inputs_after_fix, outputs_after_fix, 1000, num_epochs, learning_rate,lambda);

    // printf("Optimal parameters found: a = %.2f, b = %.2f, c = %.2f, d = %.2f,function[%d]= %p\n", neural->best_a, neural->best_b, neural->best_c, neural->best_d,0,neural->fx);


    // generate_data(sample_inputs, sample_outputs, 10);
    // for(int i=0;i<10;i++){       
    //     predicted_out[i] = forward_neural(neural, sample_inputs[i]);
    //     printf("Sample %d, x = %.3f, result =%.3f, Predicted = %.3f\n", i, sample_inputs[i], sample_outputs[i],predicted_out[i]);
    // }
    
    // compute_loss(neural, sample_inputs, sample_outputs, 10);


    free(neural->coef);
    free(neural->grid);
    free(neural);
    
    return 0;
}
