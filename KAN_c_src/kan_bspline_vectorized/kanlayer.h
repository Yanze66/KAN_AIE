#ifndef KANLAYER_H
#define KANLAYER_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "spline2fun_adam.h"


typedef struct {
    int in_dim;      // 输入维度
    int out_dim;     // 输出维度
    Neural **neurons; // 神经元数组
    double *node_entropy; // 存储节点的输出
} KANLayer;

void xavier_initialization(KANLayer* KANLayer, Neural* neural) ;

KANLayer* init_kan_layer(int layer_num, int in_dim, int out_dim, int num, int k) ;

void forward_kan_layer(KANLayer *layer, double *inputs, double *outputs) ;

void backward_kan_layer(KANLayer *layer, double *inputs, double *errors, double learning_rate , double lambda) ;
void prune_kanlayer(KANLayer *layer, double threshold) ;

#endif