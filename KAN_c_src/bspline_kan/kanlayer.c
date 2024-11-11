#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "spline2fun_adam.h"
#include "kanlayer.h"

//  define the entropy function
// double calculate_entropy(double *values, int length) {
//     double sum = 0.0;
//     double entropy = 0.0;

//     // 计算总和
//     for (int i = 0; i < length; i++) {
//         sum += values[i];
//     }

//     // 计算熵
//     for (int i = 0; i < length; i++) {
//         if (values[i] > 0) {
//             double probability = values[i] / sum;
//             entropy -= probability * log(probability);
//         }
//     }
//     return entropy;
// }

double calculate_variance(double *values, int length) {
    double sum = 0.0;
    double mean = 0.0;
    double variance = 0.0;

    // 计算均值
    for (int i = 0; i < length; i++) {
        sum += values[i];
    }
    mean = sum / length;

    // 计算方差
    for (int i = 0; i < length; i++) {
        variance += (values[i] - mean) * (values[i] - mean);
    }
    return variance / length;
}

// 更新神经元输出历史并计算熵
double update_and_calculate_entropy(Neural *neural, double current_output) {
    // 更新输出历史
    neural->output_history[neural->history_index] = current_output;
    neural->history_index = (neural->history_index + 1) % 20;

    // 计算并返回熵
    return calculate_variance(neural->output_history, 20);
}

// Xavier 初始化函数
void xavier_initialization(KANLayer* KANLayer, Neural* neural) {
    // 计算 Xavier 系数
    double scale = sqrt(6.0 / (KANLayer->in_dim + KANLayer->out_dim));

    // 使用均匀分布 -scale 到 +scale 之间初始化 wb
    neural->wb = ((double)rand() / RAND_MAX) * 2 * scale - scale;
}


KANLayer* init_kan_layer(int layer_num, int in_dim, int out_dim, int num, int k) {
    KANLayer *layer = (KANLayer*)malloc(sizeof(KANLayer));
    layer->in_dim = in_dim;
    layer->out_dim = out_dim;
    
    layer->neurons = (Neural**)malloc(in_dim * out_dim * sizeof(Neural*));
    layer->node_entropy = (double*)calloc( in_dim * out_dim, sizeof(double));  //sotre 20 nodes' entropy

    for(int i = 0; i < in_dim; i++) {
        for (int j = 0; j < out_dim; j++) {
            
            layer->neurons[i*out_dim+j] = init_neural(layer_num, i, j, num, k);
            xavier_initialization(layer, layer->neurons[i*out_dim+j]); // Xavier initialization to wb
        }
    }

    return layer;
}

void forward_kan_layer(KANLayer *layer, double *inputs, double *outputs) {
    double output[layer->out_dim * layer->out_dim];
    int num = 20; // 用于存储历史值的数量
    for (int j = 0; j < layer->out_dim; j++) {
        double output_sum = 0.0;
        for (int i = 0; i < layer->in_dim; i++) {
            // 以每个输入节点激发多个输出
            if(layer->neurons[j*layer->in_dim+i]->active){
                //printf (" |neural[%d]'s input = %.3f\n",j*layer->out_dim+i,inputs[i]);
                double node_output = forward_neural(layer->neurons[j*layer->in_dim+i], inputs[i]);
                
                layer->node_entropy[j*layer->in_dim+i] = update_and_calculate_entropy(layer->neurons[j*layer->in_dim+i], node_output);
                output_sum += node_output;
            }
     
        }
        outputs[j] = output_sum;
        //printf(" |result: %f\n",outputs[j]);
    }

}

void backward_kan_layer(KANLayer *layer, double *inputs, double *errors, double learning_rate, double lambda) {
    for (int j = 0; j < layer->out_dim; j++) { 
        for (int i = 0; i < layer->in_dim; i++) {         
            // 更新节点
            backward_neural(layer->neurons[j*layer->in_dim+i], inputs[i], errors[j], learning_rate, lambda);
        }
    }
}

// void train_single_node(Neural *neuron, double *inputs, double *outputs, int num_samples, double learning_rate, double lambda) {
//        for (int i = 0; i < num_samples; i++) {
//            double x = inputs[i];
//            double true_value = outputs[i];
//            double predicted_value = forward_neural(neuron, x);
//            double error = predicted_value - true_value;
//            backward_neural(neuron, x, error, learning_rate, lambda);
//        }
//    }




void prune_kanlayer(KANLayer *layer, double threshold) {
    for (int i = 0; i < layer->in_dim; i++) {
        for (int j = 0; j < layer->out_dim; j++) {
            Neural *neural = layer->neurons[i*layer->out_dim+j];
            if (neural->active && layer->node_entropy[i*layer->out_dim+j] < threshold) {
                neural->active = 0;
            }
        }
    }
}

