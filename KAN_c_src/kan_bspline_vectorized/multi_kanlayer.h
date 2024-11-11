#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "spline2fun_adam.h"
#include "kanlayer.h"

typedef struct {
    KANLayer **layers; // 指向层的指针数组
    int num_layers;    // 层的数量
} KANNetwork;


KANNetwork* init_kan_network(int num_layers, int *layer_dims, int num, int k) ;
void forward_kan_network(KANNetwork *network, double *inputs, double *final_output) ;

void backward_kan_network(KANNetwork *network,  double output_errors, double learning_rate, double lambda) ;

void prune_kannetwork(KANNetwork *network, double threshold) ;