#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#define  M_PI 3.1415926

// define the structure of a single neuron
typedef struct {
    double wb;   // weight of the bias
    double ws;   // weight of the spline
    double *coef; // Bspline coefficients
    double *grid; // Bspline grid
    int active;   // neuron activation status
    
} Neural;

// define the structure of a KAN layer, including multiple neurons
typedef struct {
    int in_dim, out_dim, num, k; //  input dimension, output dimension, number of splines, order of splines
    Neural **neurons; // neurons in the layer
} KANLayer;

//   define the structure of a KAN network, including multiple layers
typedef struct {
    int num_layers;
    KANLayer **layers;
} KANNetwork;

//  define the function to calculate the B-spline basis function
void B_batch(double x, double *grid, int n, int k, double *result) {
    // initialize the result to zero, totally n+k-1 basis functions
    for (int i = 0; i < n + k - 1; i++) {
        result[i] = 0.0;
    }

    //  calculate the B-spline basis function recursively
    if (k == 0) {
        for (int g = 0; g < n - 1; g++) {
            if (x >= grid[g] && x < grid[g + 1]) {
                result[g] = 1.0;
            }
        }
    } else {
        double B_km1[n + k - 2]; // bspline basis function of k-1
        B_batch(x, grid, n, k - 1, B_km1);

        for (int g = 0; g < n + k - 2; g++) {
            double alpha_denominator = grid[g + k] - grid[g];
            double beta_denominator = grid[g + k + 1] - grid[g + 1];
            double alpha = (alpha_denominator == 0) ? 0 : (x - grid[g]) / alpha_denominator;
            double beta = (beta_denominator == 0) ? 0 : (grid[g + k + 1] - x) / beta_denominator;
            result[g] = alpha * B_km1[g] + beta * B_km1[g + 1];
        }
    }
}

//  define the function to extend the grid by k points on both sides
void extend_grid(double *grid, int grid_points, int k, double *result) {
    double h = (grid[grid_points - 1] - grid[0]) / (grid_points - 1);
    for (int j = 0; j < grid_points; j++) {
        result[j + k] = grid[j];
    }
    for (int j = 0; j < k; j++) {
        result[j] = grid[0] - (k - j) * h;
        result[grid_points + k + j] = grid[grid_points - 1] + (j + 1) * h;
    }
}

double relu(double x) {
    return x > 0 ? x : 0;
}

//  define the function to initialize a single neuron
Neural* init_neural( int num, int k, double noise_scale) {
    Neural* neural = (Neural*)malloc(sizeof(Neural));
    
    //double limit = sqrt(6.0 / (in_dim + out_dim)); // Xavier initialization
    neural->wb = ((double)rand() / RAND_MAX) * noise_scale; // He initialization for wb
    
    neural->ws = ((double)rand() / RAND_MAX) * noise_scale; // He initialization for ws

    neural->coef = (double*)malloc(num * sizeof(double));
    neural->grid = (double*)malloc((num + 2 * k) * sizeof(double));
    
    for (int n = 0; n < num; n++) {
        neural->coef[n] = 1.0; // Approximate initialization for uniform output
        //neural->coef[n] = ((double)rand() / RAND_MAX - 0.5) * 0.01; // Initialize close to zero
    }
    //  initialize the grid
    double step = 1.0 / (num + 2 * k); //  uniform grid
    for (int i = 0; i < num + 2 * k; i++) {
        neural->grid[i] = i * step;
    }
    neural->active = 1; //  active neuron

    return neural;
}

//  define the function to free the memory of a single neuron
void free_neural(Neural *neural) {
    free(neural->coef);
    free(neural->grid);
    free(neural);
}
//  define the function to free the memory of a KAN layer
void free_kanlayer(KANLayer *layer) {
    for (int i = 0; i < layer->in_dim; i++) {
        for (int j = 0; j < layer->out_dim; j++) {
              free_neural(layer->neurons[i * layer->out_dim + j]);
        }
    }
    free(layer->neurons);
    free(layer);
}

//  define the function to free the memory of a KAN network
void free_kan_network(KANNetwork *network) {
    for (int l = 0; l < network->num_layers - 1; l++) {
        free_kanlayer(network->layers[l]);
    }
    free(network->layers);
    free(network);
}

//  define the function to prune the KAN layer
void prune_kanlayer(KANLayer *layer, double threshold) {
    for (int i = 0; i < layer->in_dim; i++) {
        for (int j = 0; j < layer->out_dim; j++) {
            Neural *neural = layer->neurons[i * layer->out_dim + j];
            if (fabs(neural->wb) + fabs(neural->ws) < threshold) {
                neural->active = 0;
            }
        }
    }
}

//  define the function to initialize a KAN layer, including multiple neurons, input dimension, output dimension, number of splines, order of splines, noise scale
KANLayer* init_kanlayer(int in_dim, int out_dim, int num, int k, double noise_scale) {
    KANLayer *layer = (KANLayer*)malloc(sizeof(KANLayer));
     
    layer->in_dim = in_dim;
    layer->out_dim = out_dim;
    layer->num = num;
    layer->k = k;
    layer->neurons = (Neural**)malloc(in_dim * out_dim * sizeof(Neural*));
     
    for (int i = 0; i < in_dim; i++) {
        for (int j = 0; j < out_dim; j++) {
            layer->neurons[i * out_dim + j] = init_neural(num, k, noise_scale);
        }
    }
    return layer;
}

//  define the function to calculate the forward output of a single neuron, input is the input data, num is the number of splines, k is the order of splines
double forward_neural(Neural *neural, double input,int num, int k) {
    double sx_result[num + k];
    B_batch(input, neural->grid, num, k, sx_result);

    double spline_value = 0.0;
    for (int n = 0; n < num; n++) {
        spline_value += neural->coef[n] * sx_result[n];
    }

    return neural->active * (relu(input) * neural->wb + spline_value * neural->ws);
}

//  define the function to calculate the forward output of a KAN layer,input is the input data,output is the output data
void forward_kanlayer(KANLayer *layer, double *input, double *output) {
    for (int j = 0; j < layer->out_dim; j++) {
        output[j] = 0.0;
        for (int i = 0; i < layer->in_dim; i++) {
            Neural *neural = layer->neurons[i * layer->out_dim + j];
            output[j] += forward_neural(neural, input[i],layer->num,layer->k);
        }
    }
}




//  define the function to initialize a KAN network, including multiple layers, layer sizes, number of layers, number of splines, order of splines, noise scale
KANNetwork* init_kan_network(int *layer_sizes, int num_layers, int num, int k, double noise_scale) {
    KANNetwork *network = (KANNetwork*)malloc(sizeof(KANNetwork));
    
    network->num_layers = num_layers;
    network->layers = (KANLayer**)malloc((num_layers - 1) * sizeof(KANLayer*));
     

    for (int l = 0; l < num_layers - 1; l++) {
        network->layers[l] = init_kanlayer(layer_sizes[l], layer_sizes[l + 1], num, k, noise_scale);
        
    }

    return network;
}

//  
//  define the function to calculate the backward output of a single neuron, input is the input data, num is the number of splines, k is the order of splines
// void backward_single_sample(KANNetwork *network, double *outputs, double *expected_outputs, int out_dim, int num_samples, double learning_rate) {
//     //double error_sum = 0.0;
//     double total_error_sum = 0.0;
//     for(int n=0;n<num_samples;n++){
//         for (int j = 0; j < out_dim; j++) {
//             double error = outputs[n * out_dim + j] - expected_outputs[n * out_dim + j];       
//             total_error_sum += error * error;        
//         }
//     }
    

//     //  calculate the RMSE
//     double rmse = sqrt(total_error_sum / (out_dim * num_samples));
//     if (rmse == 0) { 
//         fprintf(stderr, "Warning: RMSE is zero, potential division by zero.\n");
//         return; //  avoid division by zero
//     }

//      // update weights and B-spline coefficients through backpropagation
//     for (int n = 0; n < num_samples; n++) {
//         double loss_derivative[out_dim];
//         for (int j = 0; j < out_dim; j++) {
//             //  calculate the loss derivative
//             loss_derivative[j] = (outputs[n * out_dim + j] - expected_outputs[n * out_dim + j]) / (rmse * out_dim);
//         }

//         //  update weights and B-spline coefficients through backpropagation
//         for (int l = network->num_layers - 2; l >= 0; l--) {
//             KANLayer *layer = network->layers[l];

//             for (int i = 0; i < layer->in_dim; i++) {
//                 for (int j = 0; j < layer->out_dim; j++) {
//                     Neural *neural = layer->neurons[i * layer->out_dim + j];

//                     if (!neural->active) continue;

//                     double activation_derivative = outputs[n * out_dim + j] > 0 ? 1 : 0; // ReLU derivative

//                     //  calculate the gradient of the bias and the spline
//                     double gradient_wb = activation_derivative * loss_derivative[j];
//                     double gradient_ws = gradient_wb;

//                     //  update the weights 
//                     neural->wb -= learning_rate * gradient_wb;
//                     neural->ws -= learning_rate * gradient_ws;

//                     //  update the B-spline coefficients
//                     double sx_result[layer->num + layer->k];
//                     B_batch(outputs[n * out_dim + j], neural->grid, layer->num, layer->k, sx_result);

//                     for (int c = 0; c < layer->num; c++) {
//                         double gradient_coef = loss_derivative[j] * neural->ws * sx_result[c];
//                         neural->coef[c] -= learning_rate * gradient_coef;
//                     }
//                 }
//             }

//             prune_kanlayer(layer, 1e-3); //  prune the KAN layer
//         }
//     }
// }


void backward_single_sample(KANNetwork *network, double *outputs, double *expected_outputs, int out_dim, int num_samples, double learning_rate) {
    //double error_sum = 0.0;
    double total_error_sum = 0.0;
    for(int n=0;n<num_samples;n++){
        for (int j = 0; j < out_dim; j++) {
            double error = outputs[n * out_dim + j] - expected_outputs[n * out_dim + j];       
            total_error_sum += error * error;        
        }
    }
    

    //  calculate the RMSE
    double rmse = sqrt(total_error_sum / (out_dim * num_samples));
    if (rmse == 0) { 
        fprintf(stderr, "Warning: RMSE is zero, potential division by zero.\n");
        return; //  avoid division by zero
    }

     // update weights and B-spline coefficients through backpropagation
    for (int n = 0; n < num_samples; n++) {
        double loss_derivative[out_dim];
        for (int j = 0; j < out_dim; j++) {
            //  calculate the loss derivative
            loss_derivative[j] = (outputs[n * out_dim + j] - expected_outputs[n * out_dim + j]) / (rmse * out_dim); ;
        }

        //  update weights and B-spline coefficients through backpropagation
        for (int l = network->num_layers - 2; l >= 0; l--) {
            KANLayer *layer = network->layers[l];

            for (int i = 0; i < layer->in_dim; i++) {
                for (int j = 0; j < layer->out_dim; j++) {
                    Neural *neural = layer->neurons[i * layer->out_dim + j];

                    if (!neural->active) continue;

                    double activation_derivative = outputs[n * out_dim + j] > 0 ? 1 : 0; // ReLU derivative

                    //  calculate the gradient of the bias and the spline
                    double gradient_wb = activation_derivative * loss_derivative[j];
                    double gradient_ws = gradient_wb;

                    //  update the weights 
                    neural->wb -= learning_rate * gradient_wb;
                    neural->ws -= learning_rate * gradient_ws;

                    //  update the B-spline coefficients
                    double sx_result[layer->num + layer->k];
                    B_batch(outputs[n * out_dim + j], neural->grid, layer->num, layer->k, sx_result);

                    for (int c = 0; c < layer->num; c++) {
                        double gradient_coef = loss_derivative[j] * neural->ws * sx_result[c];
                        neural->coef[c] -= learning_rate * gradient_coef;
                    }
                }
            }

            prune_kanlayer(layer, 1e-3); //  prune the KAN layer
        }
    }
}


void forward_kan_network(KANNetwork *network, double *inputs, double *outputs, int num_samples) {
    int in_dim = network->layers[0]->in_dim;
    int out_dim = network->layers[network->num_layers - 2]->out_dim;

    //  calculate the forward output of the KAN network
    for (int n = 0; n < num_samples; n++) {
        int current_in_dim = in_dim;

        double *current_input = (double*)malloc(current_in_dim * sizeof(double));
        if (!current_input) {
            fprintf(stderr, "Error: Memory allocation for current_input failed.\n");
            exit(EXIT_FAILURE);
        }
        //  initialize the input data
        for (int i = 0; i < current_in_dim; i++) {
            current_input[i] = inputs[n * current_in_dim + i];
        }

        
        
        for (int l = 0; l < network->num_layers - 1; l++) {
            KANLayer *layer = network->layers[l];
       
            //  adjust the input dimension for the next layer
            double *current_output = (double*)malloc(layer->out_dim * sizeof(double));
            
            forward_kanlayer(layer, current_input, current_output);

            if (l < network->num_layers - 2) {
                current_input = (double*)realloc(current_input, network->layers[l + 1]->in_dim * sizeof(double));
                if (!current_input) {
                    fprintf(stderr, "Error: Memory allocation for current_input failed.\n");
                    free(current_output);
                    exit(EXIT_FAILURE);
                }
                for (int j = 0; j < network->layers[l + 1]->in_dim; j++) {
                    current_input[j] = current_output[j];
                }
            }

            if(l == network->num_layers - 2){
                //  store the output data of the KAN network
                for (int j = 0; j < out_dim; j++) {
                    outputs[n * out_dim + j] = current_output[j];
                }
            }
                    free(current_output);

        }

    
        //  free the memory of the input data
        
        free(current_input);
        
    }
}




//  generate the training data
void generate_data(double *inputs, double *outputs, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        double x = ((double)rand() / RAND_MAX) * 2 - 1; // -1 to 1
        double y = ((double)rand() / RAND_MAX) * 2 - 1;

        inputs[2 * i] = x;
        inputs[2 * i + 1] = y;

        outputs[i] = exp(sin(M_PI * x) + y * y); // test f(x, y) = exp(sin(pi*x) + y^2)
        
    }
}
//  the main function
int main() {
    srand(time(NULL)); //  set the random seed

    int layer_sizes[] = {2, 10, 1}; //  input layer, hidden layer, output layer
    int num_layers = sizeof(layer_sizes) / sizeof(layer_sizes[0]);
    int num_samples = 100; //  number of samples
    int num_epochs = 20;
    double learning_rate = 0.001;
    int spline_nodes = 10;
    int spline_degree = 3;
    double noise_scale = 0.01; //  noise scale

    //  initialize the KAN network
    KANNetwork *network = init_kan_network(layer_sizes, num_layers, spline_nodes, spline_degree, noise_scale);

    //  allocate memory for the input data, output data, and predicted output data
    double *inputs = (double*)malloc(2 * num_samples * sizeof(double)); //  input data x, y
    double *outputs = (double*)malloc(num_samples * sizeof(double));
    double *predicted_outputs = (double*)malloc(num_samples * sizeof(double));

    //  train the KAN network
    for (int epoch = 0; epoch < num_epochs; epoch++) {
        generate_data(inputs, outputs, num_samples);
        //  calculate the forward output of the KAN network
        forward_kan_network(network, inputs, predicted_outputs, num_samples);
        //  calculate the backward output of the KAN network
        backward_single_sample(network, predicted_outputs, outputs, layer_sizes[num_layers - 1], num_samples, learning_rate);
        
        //  calculate the RMSE loss
        double total_loss = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double error = predicted_outputs[i] - outputs[i];
            total_loss += error * error;
        }
        total_loss = sqrt(total_loss / num_samples);
        printf("Epoch %d: RMSE Loss = %f\n", epoch + 1, total_loss);

        
    }

    generate_data(inputs, outputs, num_samples);
    forward_kan_network(network, inputs, predicted_outputs, num_samples);
    double total_loss = 0.0;
    for (int i = 0; i < num_samples; i++) {
            double error = predicted_outputs[i] - outputs[i];
            total_loss += error * error;
        }
    total_loss = sqrt(total_loss / num_samples);
    printf(" Final RMSE Loss = %f\n", total_loss);


    //  free the memory of the input data, output data, and predicted output data
    free(inputs);
    free(outputs);
    free(predicted_outputs);
    free_kan_network(network);

    return 0;
}

