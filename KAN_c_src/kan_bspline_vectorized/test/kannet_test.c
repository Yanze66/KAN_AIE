#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../utils.h"
#include "../spline2fun_adam.h"
#include "../kanlayer.h"
#include "../multi_kanlayer.h"
void generate_3_1(double *inputs, double *outputs, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        double x = ((double)rand() / RAND_MAX) * 8-4; // -5 to 5
        double y = ((double)rand() / RAND_MAX) * 8-4;
        double z = ((double)rand() / RAND_MAX) * 8-4;

        inputs[3 * i] = x;
        inputs[3 * i + 1] = y;
        inputs[3 * i + 2] = z;

        outputs[i] = exp( 0.2*square(y) + 0.8*sin(3.1415926*z)-1); // test = exp(0.2*y + 0.8*sin(3.1415926*z)-1)
    }
}

void fit_network(KANNetwork* network){
    int fix_samples = 50;
    double inputs_for_fix[3 * fix_samples];
    double outputs_for_fix[fix_samples];
    double inputs_1_0[fix_samples], outputs_1_0[fix_samples];
    double inputs_1_1[fix_samples], outputs_1_1[fix_samples];
    double inputs_1_2[fix_samples], outputs_1_2[fix_samples];
    double inputs_1_3[fix_samples], outputs_1_3[fix_samples];
    double inputs_1_4[fix_samples], outputs_1_4[fix_samples];

    double inputs_0_0[fix_samples], outputs_0_0[fix_samples];
    double inputs_0_1[fix_samples], outputs_0_1[fix_samples];
    double inputs_0_2[fix_samples], outputs_0_2[fix_samples];
    double inputs_0_3[fix_samples], outputs_0_3[fix_samples];
    double inputs_0_4[fix_samples], outputs_0_4[fix_samples];
    double inputs_0_5[fix_samples], outputs_0_5[fix_samples];
    double inputs_0_6[fix_samples], outputs_0_6[fix_samples];
    double inputs_0_7[fix_samples], outputs_0_7[fix_samples];
    double inputs_0_8[fix_samples], outputs_0_8[fix_samples];
    double inputs_0_9[fix_samples], outputs_0_9[fix_samples];
    double inputs_0_10[fix_samples], outputs_0_10[fix_samples];
    double inputs_0_11[fix_samples], outputs_0_11[fix_samples];
    double inputs_0_12[fix_samples], outputs_0_12[fix_samples];
    double inputs_0_13[fix_samples], outputs_0_13[fix_samples];
    double inputs_0_14[fix_samples], outputs_0_14[fix_samples];
    // 生成数据
    generate_3_1(inputs_for_fix, outputs_for_fix, fix_samples);

    // 前向传播与记录输出
    for (int i = 0; i < fix_samples; i++) {
        forward_kan_network(network, &inputs_for_fix[i * 3], &outputs_for_fix[i]);

        // 记录每个节点的输入和输出
        inputs_1_0[i] = network->layers[1]->neurons[0]->in_value;
        inputs_1_1[i] = network->layers[1]->neurons[1]->in_value;
        inputs_1_2[i] = network->layers[1]->neurons[2]->in_value;
        inputs_1_3[i] = network->layers[1]->neurons[3]->in_value;
        inputs_1_4[i] = network->layers[1]->neurons[4]->in_value;
        outputs_1_0[i] = network->layers[1]->neurons[0]->out_value;
        outputs_1_1[i] = network->layers[1]->neurons[1]->out_value;
        outputs_1_2[i] = network->layers[1]->neurons[2]->out_value;
        outputs_1_3[i] = network->layers[1]->neurons[3]->out_value;
        outputs_1_4[i] = network->layers[1]->neurons[4]->out_value;

        inputs_0_0[i] = inputs_for_fix[i * 3];
        inputs_0_1[i] = inputs_for_fix[i * 3 + 1];
        inputs_0_2[i] = inputs_for_fix[i * 3 + 2];
        outputs_0_0[i] = network->layers[0]->neurons[0]->out_value;
        outputs_0_1[i] = network->layers[0]->neurons[1]->out_value;
        outputs_0_2[i] = network->layers[0]->neurons[2]->out_value;

        inputs_0_3[i] = inputs_for_fix[i * 3];
        inputs_0_4[i] = inputs_for_fix[i * 3 + 1];
        inputs_0_5[i] = inputs_for_fix[i * 3 + 2];
        outputs_0_3[i] = network->layers[0]->neurons[3]->out_value;
        outputs_0_4[i] = network->layers[0]->neurons[4]->out_value;
        outputs_0_5[i] = network->layers[0]->neurons[5]->out_value;

        inputs_0_6[i] = inputs_for_fix[i * 3];
        inputs_0_7[i] = inputs_for_fix[i * 3 + 1];
        inputs_0_8[i] = inputs_for_fix[i * 3 + 2];
        outputs_0_6[i] = network->layers[0]->neurons[6]->out_value;
        outputs_0_7[i] = network->layers[0]->neurons[7]->out_value;
        outputs_0_8[i] = network->layers[0]->neurons[8]->out_value;

        inputs_0_9[i] = inputs_for_fix[i * 3];
        inputs_0_10[i] = inputs_for_fix[i * 3 + 1];
        inputs_0_11[i] = inputs_for_fix[i * 3 + 2];
        outputs_0_9[i] = network->layers[0]->neurons[9]->out_value;
        outputs_0_10[i] = network->layers[0]->neurons[10]->out_value;
        outputs_0_11[i] = network->layers[0]->neurons[11]->out_value;

        inputs_0_12[i] = inputs_for_fix[i * 3];
        inputs_0_13[i] = inputs_for_fix[i * 3 + 1];
        inputs_0_14[i] = inputs_for_fix[i * 3 + 2];
        outputs_0_12[i] = network->layers[0]->neurons[12]->out_value;
        outputs_0_13[i] = network->layers[0]->neurons[13]->out_value;
        outputs_0_14[i] = network->layers[0]->neurons[14]->out_value;

    }
    
    
    fix_fun(network->layers[1]->neurons[0], inputs_1_0, outputs_1_0, fix_samples);
    fix_fun(network->layers[1]->neurons[1], inputs_1_1, outputs_1_1, fix_samples);
    fix_fun(network->layers[1]->neurons[2], inputs_1_2, outputs_1_2, fix_samples);
    fix_fun(network->layers[1]->neurons[3], inputs_1_3, outputs_1_3, fix_samples);
    fix_fun(network->layers[1]->neurons[4], inputs_1_4, outputs_1_4, fix_samples);

    fix_fun(network->layers[0]->neurons[0], inputs_0_0, outputs_0_0, fix_samples);
    fix_fun(network->layers[0]->neurons[1], inputs_0_1, outputs_0_1, fix_samples);
    fix_fun(network->layers[0]->neurons[2], inputs_0_2, outputs_0_2, fix_samples);
    fix_fun(network->layers[0]->neurons[3], inputs_0_3, outputs_0_3, fix_samples);
    fix_fun(network->layers[0]->neurons[4], inputs_0_4, outputs_0_4, fix_samples);
    fix_fun(network->layers[0]->neurons[5], inputs_0_5, outputs_0_5, fix_samples);
    fix_fun(network->layers[0]->neurons[6], inputs_0_6, outputs_0_6, fix_samples);
    fix_fun(network->layers[0]->neurons[7], inputs_0_7, outputs_0_7, fix_samples);
    fix_fun(network->layers[0]->neurons[8], inputs_0_8, outputs_0_8, fix_samples);
    fix_fun(network->layers[0]->neurons[9], inputs_0_9, outputs_0_9, fix_samples);
    fix_fun(network->layers[0]->neurons[10], inputs_0_10, outputs_0_10, fix_samples);
    fix_fun(network->layers[0]->neurons[11], inputs_0_11, outputs_0_11, fix_samples);
    fix_fun(network->layers[0]->neurons[12], inputs_0_12, outputs_0_12, fix_samples);
    fix_fun(network->layers[0]->neurons[13], inputs_0_13, outputs_0_13, fix_samples);
    fix_fun(network->layers[0]->neurons[14], inputs_0_14, outputs_0_14, fix_samples);


}

int layer_dimention_update(KANLayer* layers){
    int active_num = 0;
    for(int i = 0; i < layers->in_dim; i++){
        layers->neurons[i]->active = 1;
        active_num++;
    }
    return active_num;
}

int main() {
    srand(time(NULL)); //  set the random seed
    int num_layers = 2;
    int layer_dims[] = {3, 5, 1};
    int num = 20;
    int k = 3;
    int num_samples = 10000;
    int num_epochs = 10;
    double learning_rate = 0.001;
    double lambda = 0.001; // regularization parameter

    // 初始化两层 KAN 网络
    KANNetwork* network = init_kan_network(num_layers, layer_dims, num, k);

    // 设置数据集
    double inputs[3*num_samples];
    double outputs[num_samples];
    generate_3_1(inputs, outputs, num_samples);
    
   
    // 训练网络
    double final_output[num_samples];
    //train and test and prune the network. its a loop until find the best branch.
    for(int prune_num=0;prune_num <20; prune_num++){
        //initialize the active neural node, otherwise the network will be pruned to zero.
        // if(prune_num  == 10){
        //     for(int i = 0; i < num_layers; i++){
        //         if(i==0){
        //             network->layers[i]->in_dim = 3;                   
        //         }else{
        //             network->layers[i]->in_dim = layer_dimention_update(network->layers[i]);
        //         }
        //         for(int j = 0; j < network->layers[i]->out_dim; j++){
        //             for(int h = 0; h < network->layers[i]->in_dim; h++){
        //                 if(network->layers[i]->neurons[j*network->layers[i]->in_dim+h]->active){
        //                     network->layers[i]->neurons[j*network->layers[i]->in_dim+h] = init_neural(i, j*network->layers[i]->in_dim+h, j, num, k);
        //                     xavier_initialization(network->layers[i], network->layers[i]->neurons[j*network->layers[i]->in_dim+h]); // Xavier initialization to wb
        //                     printf("initialize the neural[%d] in layer[%d]\n",j*network->layers[i]->in_dim+h,i);
        //                 }
        //             }   
        //         }
            
        //     }
        // }
        
        printf("\ntrain the network in %d rounds\n",prune_num);
        for (int epoch = 0; epoch < num_epochs; epoch++) {
            for (int i = 0; i < num_samples; i++) {
                double new_learning_rate=0;
                forward_kan_network(network, &inputs[i * 3], &final_output[i]);
                double error = final_output[i] - outputs[i];
                // if(i%1000 == 0){
                //     printf("layer_inputs = %.3f, %.3f, %.3f, final_output = %.3f, outputs = %.3f\n", network->layers[1]->neurons[0]->in_value, network->layers[1]->neurons[1]->in_value, network->layers[1]->neurons[2]->in_value, final_output[i], outputs[i]);
                //     //printf("error = %f, final_output = %f, outputs = %f\n", error, final_output[i], outputs[i]);
                // }
                
                if(network->layers[1]->neurons[0]->active + network->layers[1]->neurons[1]->active + network->layers[1]->neurons[2]->active + network->layers[1]->neurons[3]->active + network->layers[1]->neurons[4]->active == 0){
                    printf("All nodes in layer[1] are pruned, break\n");
                    break;
                }
                
                new_learning_rate = learning_rate/(1+0.1*prune_num); ;
                backward_kan_network(network, error, new_learning_rate, lambda);
            }
            double rmse = calculate_rmse(outputs, final_output, num_samples);
            printf("Epoch %d: RMSE Loss = %f\n", epoch + 1, rmse);
        }


        // 测试网络
        printf("Testing network after training\n"); 
        generate_3_1(inputs, outputs, 50);
        for (int i = 0; i < 50; i++) {
            forward_kan_network(network, &inputs[i * 3], &final_output[i]);
            // if(i%10 == 0){
            //         printf("layer_inputs = %.3f, %.3f, %.3f, %.3f, %.3f, layer_outputs = %.3f, %.3f, %.3f, %.3f, %.3f, final_output = %.3f, outputs = %.3f\n", network->layers[1]->neurons[0]->in_value, network->layers[1]->neurons[1]->in_value, network->layers[1]->neurons[2]->in_value, network->layers[1]->neurons[3]->in_value, network->layers[1]->neurons[4]->in_value, network->layers[1]->neurons[0]->out_value, network->layers[1]->neurons[1]->out_value, network->layers[1]->neurons[2]->out_value, network->layers[1]->neurons[3]->out_value, network->layers[1]->neurons[4]->out_value, final_output[i], outputs[i]);
            //     }
        }
        double rmse_final = calculate_rmse(outputs, final_output, 50);
        printf("Final RMSE Loss : %f\n", rmse_final);


        // 修剪网络
        //prune_kannetwork(network, 0.1);
        // prune_kanlayer(network->layers[1], 0.1);
        // if(prune_num %3 ==2){
        //     prune_kanlayer(network->layers[0], 0.01);
        // }
        // printf("node in layer[1]: neural[0]->active = %d, neural[1]->active = %d, neural[2]->active = %d,neural[3]->active = %d,neural[4]->active = %d\n", network->layers[1]->neurons[0]->active, network->layers[1]->neurons[1]->active, network->layers[1]->neurons[2]->active, network->layers[1]->neurons[3]->active, network->layers[1]->neurons[4]->active);
        // printf("node in layer[0]:%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",network->layers[0]->neurons[0]->active,network->layers[0]->neurons[1]->active,network->layers[0]->neurons[2]->active,network->layers[0]->neurons[3]->active,network->layers[0]->neurons[4]->active,network->layers[0]->neurons[5]->active,network->layers[0]->neurons[6]->active,network->layers[0]->neurons[7]->active,network->layers[0]->neurons[8]->active,network->layers[0]->neurons[9]->active,network->layers[0]->neurons[10]->active,network->layers[0]->neurons[11]->active,network->layers[0]->neurons[12]->active,network->layers[0]->neurons[13]->active,network->layers[0]->neurons[14]->active);
    }
    
    
    //fit the nerual network
    // fit_network(network);

    generate_3_1(inputs, outputs, 50);
        for (int i = 0; i < 50; i++) {
            forward_kan_network(network, &inputs[i * 3], &final_output[i]);
            if(i%10 == 0){
                    printf("layer_inputs = %.3f, %.3f, %.3f, %.3f, %.3f, layer_outputs = %.3f, %.3f, %.3f, %.3f, %.3f, final_output = %.3f, outputs = %.3f\n", network->layers[1]->neurons[0]->in_value, network->layers[1]->neurons[1]->in_value, network->layers[1]->neurons[2]->in_value, network->layers[1]->neurons[3]->in_value, network->layers[1]->neurons[4]->in_value, network->layers[1]->neurons[0]->out_value, network->layers[1]->neurons[1]->out_value, network->layers[1]->neurons[2]->out_value, network->layers[1]->neurons[3]->out_value, network->layers[1]->neurons[4]->out_value, final_output[i], outputs[i]);
                }
        }
        double rmse_final = calculate_rmse(outputs, final_output, 50);
        printf("Final RMSE Loss after pruning: %f\n", rmse_final);


   
    free(network);

    return 0;
}