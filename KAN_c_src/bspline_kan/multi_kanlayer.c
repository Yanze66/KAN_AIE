#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "spline2fun_adam.h"
#include "kanlayer.h"
#include "multi_kanlayer.h"

// typedef struct {
//     KANLayer **layers; // 指向层的指针数组
//     int num_layers;    // 层的数量
// } KANNetwork;


KANNetwork* init_kan_network(int num_layers, int *layer_dims, int num, int k) {
    KANNetwork *network = (KANNetwork*)malloc(sizeof(KANNetwork));
    network->num_layers = num_layers;
    network->layers = (KANLayer**)malloc(num_layers * sizeof(KANLayer*));

    for (int i = 0; i < num_layers; i++) {
        network->layers[i] = init_kan_layer(i, layer_dims[i], layer_dims[i + 1], num, k);
    }

    return network;
}

void forward_kan_network(KANNetwork *network, double *inputs, double *final_output) {
    double *current_inputs = inputs;
    double *layer_outputs = NULL;

    for (int i = 0; i < network->num_layers; i++) {
        int layer_output_dim = network->layers[i]->out_dim;
        layer_outputs = (double *)malloc(layer_output_dim * sizeof(double));
        //printf("Layer =%d, input = %.3f, %.3f, %.3f,layer_output_dim = %d,layer_outputs = %.3f, %.3f, %.3f\n", i, current_inputs[0], current_inputs[1], current_inputs[2], layer_output_dim, layer_outputs[0], layer_outputs[1], layer_outputs[2]);
        forward_kan_layer(network->layers[i], current_inputs, layer_outputs);
        if(i == 0){
            //printf("Layer =%d, input = %.3f, %.3f, %.3f, output = %.3f, %.3f, %.3f,%.3f, %.3f\n", i, current_inputs[0], current_inputs[1], current_inputs[2], layer_outputs[0], layer_outputs[1], layer_outputs[2],layer_outputs[3], layer_outputs[4]);
        }
        if(i ==1){
            // printf("Layer =%d, input = %.3f, %.3f, %.3f, %.3f, %.3f, output = %.3f\n", i, current_inputs[0], current_inputs[1], current_inputs[2],current_inputs[3],current_inputs[4], layer_outputs[0]);
        }
        current_inputs = layer_outputs; // 将这一层的输出作为下一层的输入
    }


    for (int i = 0; i < network->layers[network->num_layers - 1]->out_dim; i++) {
        final_output[i] = layer_outputs[i];
    }

    free(layer_outputs); // 释放最后一次的输出缓冲区
}

void backward_kan_network(KANNetwork *network,  double output_errors, double learning_rate, double lambda) {
    // 
    double *errors = NULL;
    double *inputs = NULL;


    for (int i = network->num_layers - 1; i >= 0; i--) {
        inputs = (double *)malloc(network->layers[i]->in_dim * sizeof(double));
        errors = (double *)malloc(network->layers[i]->out_dim * sizeof(double)); 

        if(i == network->num_layers - 1){
            errors[0] = output_errors;
            //printf("Layer =%d, output_errors = %.3f, errors = %.3f\n", i, output_errors, errors[0]);
        }else{
            for (int h = 0; h < network->layers[i+1]->in_dim; h++) {
                errors[h] = network->layers[i+1]->neurons[h]->pre_error;
            }
            //printf("Layer =%d, output_errors = %.3f, errors = %.3f, %.3f, %.3f\n", i, output_errors, errors[0], errors[1], errors[2]);
        } 

        for(int j = 0; j < network->layers[i]->out_dim; j++){   
            for (int k = 0; k < network->layers[i]->in_dim; k++) {               
                if(network->layers[i]->neurons[k+network->layers[i]->in_dim*j]->active){
                    inputs[k] = network->layers[i]->neurons[k+network->layers[i]->in_dim*j]->in_value;
                    calculate_previous_error(network->layers[i]->neurons[k+network->layers[i]->in_dim*j], errors[j]);             
                }else{
                    inputs[k] = 0;
                    network->layers[i]->neurons[k+network->layers[i]->in_dim*j]->pre_error = 0;
                    //calculate_previous_error(network->layers[i]->neurons[k+network->layers[i]->in_dim*j], errors[j]);             

                }
                //inputs[k] = network->layers[i]->neurons[k+network->layers[i]->in_dim*j]->in_value;
                //calculate_previous_error(network->layers[i]->neurons[k+network->layers[i]->in_dim*j], errors[j]);             
            }        
            
        }
        
        
        //printf("Layer =%d, neuron = %d, input = %.3f, error = %.3f\n", i, k*network->layers[i]->out_dim+j, inputs[k], errors);
        backward_kan_layer(network->layers[i], inputs, errors, learning_rate, lambda);
        
        free(inputs);
        free(errors);
        
        
    }
}

void prune_kannetwork(KANNetwork *network, double threshold) {
    //prune_kanlayer(network->layers[network->num_layers-1], threshold);
    
    for (int i = 0; i < network->num_layers; i++) {
        prune_kanlayer(network->layers[i], threshold);
    
    }
}
