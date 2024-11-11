#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../utils.h"
#include "../spline2fun_adam.h"
#include "../kanlayer.h"
generate_3_1(double *inputs, double *outputs, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        double x = ((double)rand() / RAND_MAX) * 8-4; // -5 to 5
        double y = ((double)rand() / RAND_MAX) * 8-4;
        double z = ((double)rand() / RAND_MAX) * 8-4;

        inputs[3 * i] = x;
        inputs[3 * i + 1] = y;
        inputs[3 * i + 2] = z;

        outputs[i] = 1 - 2*square(y) + 2*sin(3.1415926*z); // test = 2*sin(3.1415926*x) - 6*y^2 + sqrt(z)
    }
}

void train_layer(KANLayer *layer, double *inputs, double *outputs, int num_samples, int num_epochs, double learning_rate, double lambda) {
    double *layer_outputs = (double *)malloc(layer->out_dim *num_samples * sizeof(double));
    double *output_errors = (double *)malloc(layer->out_dim *num_samples * sizeof(double)); // 确保正确大小

    for (int epoch = 0; epoch < num_epochs; epoch++) {
        double total_loss = 0.0;
        
        for (int i = 0; i < num_samples; i++) {
            // 前向传播
            forward_kan_layer(layer, &inputs[i*3], &layer_outputs[i]);
            double error[] = {0.0}; //only 1 dimention for test.
            error[0] = layer_outputs[i] - outputs[i];
            total_loss += error[0] * error[0];

            // 记录每个输出节点的误差
            for (int j = 0; j < layer->out_dim; j++) {
                output_errors[j] = error[j]; // 如果每个输出有其特定误差来源，这里可以调整
            }

            // 反向传播
            backward_kan_layer(layer, &inputs[i*3], error, learning_rate, lambda);
        }

        total_loss = sqrt(total_loss / num_samples);
        printf("Epoch %d: RMSE Loss = %f\n", epoch + 1, total_loss);
    }
    
    free(layer_outputs);
    free(output_errors);
}

void test(KANLayer *layer) {
    int num_samples = 50;
    double test_inputs[3 * num_samples];
    double test_outputs[num_samples];
    double predicted_outputs[num_samples];

    // 生成测试数据
    generate_3_1(test_inputs, test_outputs, num_samples);

    // 使用测试数据进行预测
    for (int i = 0; i < num_samples; i++) {
        forward_kan_layer(layer, &test_inputs[i * 3], &predicted_outputs[i]);

        // for(int j = 0; j < layer->out_dim * layer->in_dim; j++){
        //     printf(" |node_entropy: %f\n",layer->node_entropy[j]);
        // }
        printf("Sample %d, x = %.3f,y = %.3f,z = %.3f, neural[0]= %.3f, neural[1]= %.3f, neural[2]= %.3f, Predicted = %.3f, Real = %.3f\n", i, test_inputs[i * 3], test_inputs[i * 3 + 1], test_inputs[i * 3 + 2], layer->neurons[0]->out_value, layer->neurons[1]->out_value, layer->neurons[2]->out_value, predicted_outputs[i], test_outputs[i]);
    }
    // 计算并输出 RMSE
    double rmse = calculate_rmse(test_outputs, predicted_outputs, num_samples);
    printf("Final RMSE Loss: %f\n", rmse);
}


void fit_layer(KANLayer *layer) {
    int fix_samples = 50;
    double inputs_for_fix[3 * fix_samples];
    double outputs_for_fix[fix_samples];
    double inputs1[fix_samples], outputs1[fix_samples];
    double inputs2[fix_samples], outputs2[fix_samples];
    double inputs3[fix_samples], outputs3[fix_samples];

    // 生成数据
    generate_3_1(inputs_for_fix, outputs_for_fix, fix_samples);

    // 前向传播与记录输出
    for (int i = 0; i < fix_samples; i++) {
        forward_kan_layer(layer, &inputs_for_fix[i * 3], &outputs_for_fix[i]);

        // 记录每个节点的输入和输出
        inputs1[i] = inputs_for_fix[i * 3];
        inputs2[i] = inputs_for_fix[i * 3 + 1];
        inputs3[i] = inputs_for_fix[i * 3 + 2];
        outputs1[i] = layer->neurons[0]->out_value;
        outputs2[i] = layer->neurons[1]->out_value;
        outputs3[i] = layer->neurons[2]->out_value;
    }

    // 函数替换
    fix_fun(layer->neurons[0], inputs1, outputs1, fix_samples);
    fix_fun(layer->neurons[1], inputs2, outputs2, fix_samples);
    fix_fun(layer->neurons[2], inputs3, outputs3, fix_samples);
}


int main() {
    srand(time(NULL)); //  set the random seed
    int in_dim = 3;
    int out_dim = 1;
    int num = 50;
    int k = 3;
    int num_samples = 20000;
    int num_epochs = 100;
    double learning_rate = 0.001;
    double lambda = 0.001; // regularization parameter

    KANLayer *layer = init_kan_layer(0, in_dim,out_dim, num, k);

    // 设置数据集
    double inputs[3*num_samples];
    double outputs[num_samples];
    generate_3_1(inputs, outputs, num_samples);
    //standardize_data(inputs, 3*num_samples);
    //standardize_data(outputs, num_samples);

    // 定义输出层的误差存储
    double layer_outputs[out_dim*num_samples];
    double output_errors[out_dim*num_samples];

    
    // 训练网络
    train_layer(layer, inputs, outputs, num_samples, num_epochs, learning_rate, lambda);

    

    // 测试网络
    test(layer);

    prune_kanlayer(layer, 0.1);

     // 训练网络
    generate_3_1(inputs, outputs, num_samples);
    train_layer(layer, inputs, outputs, num_samples, num_epochs, learning_rate, lambda);

    

    // 测试网络
    test(layer);

    // 修正网络
    fit_layer(layer);

    

    // // 重新测试网络
    // test(layer);



    // printf("layer->neurons[0]->best_a = %.2f, layer->neurons[0]->best_b = %.2f, layer->neurons[0]->best_c = %.2f, layer->neurons[0]->best_d = %.2f\n", layer->neurons[0]->best_a, layer->neurons[0]->best_b, layer->neurons[0]->best_c, layer->neurons[0]->best_d);
    // printf("layer->neurons[1]->best_a = %.2f, layer->neurons[1]->best_b = %.2f, layer->neurons[1]->best_c = %.2f, layer->neurons[1]->best_d = %.2f\n", layer->neurons[1]->best_a, layer->neurons[1]->best_b, layer->neurons[1]->best_c, layer->neurons[1]->best_d);
    // printf("layer->neurons[2]->best_a = %.2f, layer->neurons[2]->best_b = %.2f, layer->neurons[2]->best_c = %.2f, layer->neurons[2]->best_d = %.2f\n", layer->neurons[2]->best_a, layer->neurons[2]->best_b, layer->neurons[2]->best_c, layer->neurons[2]->best_d);

    //重新train
    // train_layer(layer, inputs, outputs, num_samples, num_epochs, learning_rate, lambda);

    // 重新测试
    test(layer);

    // //free(layer->neurons);
     free(layer);

    return 0;
}