
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
double identity(double x) {
    return x;
}

double square(double x) {
    return x * x;
}

double cube(double x) {
    return x * x * x;
}


double power4(double x) {
    return pow(x, 4);
}

double power5(double x) {
    return pow(x, 5);
}


double reciprocal(double x) {
    if (x != 0)
        return 1.0 / x;
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal function.\n");
        return 100;
}

double reciprocal_square(double x) {
    if (x != 0)
        return 1.0 / (x * x);
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal_square function.\n");
        return 100;
}

double reciprocal_cube(double x) {
    if (x != 0)
        return 1.0 / (x * x * x);
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal_cube function.\n");
        return 100;
}

double reciprocal_x4(double x) {
    if (x != 0)
        return 1.0 / (x * x * x * x);
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal_cube function.\n");
        return 100;
}

double reciprocal_x5(double x) {
    if (x != 0)
        return 1.0 / (x * x * x *x *x);
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal_cube function.\n");
        return 100;
}

//y = sqrt(x)
double sqrt_func(double x) {
    if (x >= 0)
        return sqrt(x);
    else
        // fprintf(stderr, "Warning: Negative input to sqrt function.\n");
        return 100;
}

//y =sqrt(x) * sqrt(x) * sqrt(x)
double cube_root(double x) {
    if (x >= 0)
        return sqrt(x) * sqrt(x) * sqrt(x);
    else
        // fprintf(stderr, "Warning: Negative input to sqrt function.\n");
        return 100;
}

//y= 1/sqrt(x)
double reciprocal_sqrt(double x) {
    if (x > 0)
        return 1.0 / sqrt(x);
    else
        // fprintf(stderr, "Warning: Division by zero in reciprocal_sqrt function.\n");
        return 100;
}

//y = exp(x)
double exp_func(double x) {
    return exp(x);
}

//y=log(x)
double log_func(double x) {
    if (x > 0)
        return log(x);
    else
        // fprintf(stderr, "Warning: Negative input to log function.\n");
        return 100;
}
//y=abs(x)
double abs_func(double x) {
    return fabs(x);
}

//y=sin(x)
double sin_func(double x) {
    return sin(x);
}
//y=cos(x)
double cos_func(double x) {
    return cos(x);
}
//y=tan(x)  
double tan_func(double x) {
    return tan(x);
}
//y=tanh(x)
double tanh_func(double x) {
    return tanh(x);
}


//y=sgn(x)
double sgn_func(double x) {
    if (x > 0)
        return 1;
    else if (x < 0)
        return -1;
    else
        return 0;
}
//y=arcsin(x)
double arcsin_func(double x) {
    if (x >= -1 && x <= 1)
        return asin(x);
    else
        // fprintf(stderr, "Warning: Input out of range in arcsin function.\n");
        return 100;
}
//y=arccos(x)
double arccos_func(double x) {
    if (x >= -1 && x <= 1)
        return acos(x);
    else
        // fprintf(stderr, "Warning: Input out of range in arccos function.\n");
        return 100;
}
//y=arctan(x)
double arctan_func(double x) {
    return atan(x);
}
//y=arctanh(x)
double arctanh_func(double x) {
    if (x > -1 && x < 1)
        return atanh(x);
    else
        //fprintf(stderr, "Warning: Input out of range in arctanh function.\n");
        return 100;
}
//y=0
double zero_func(double x) {
    return 0;
}
//y=guauss(x)   
double gauss_func(double x) {
    return exp(-x * x);
}


// Derivative of y = x, i.e., identity function
double identity_derivative(double x) {
    return 1.0;
}

// Derivative of y = x^2
double square_derivative(double x) {
    return 2.0 * x;
}

// Derivative of y = x^3
double cube_derivative(double x) {
    return 3.0 * x * x;
}

// Derivative of y = x^4
double power4_derivative(double x) {
    return 4.0 * pow(x, 3);
}

// Derivative of y = x^5
double power5_derivative(double x) {
    return 5.0 * pow(x, 4);
}

// Derivative of y = 1/x
double reciprocal_derivative(double x) {
    if (x != 0) {
        return -1.0 / (x * x);
    } else {
        // Handle possible division by zero
        return 0.0;
    }
}

// Derivative of y = 1/x^2
double reciprocal_square_derivative(double x) {
    if (x != 0) {
        return -2.0 / (x * x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = 1/x^3
double reciprocal_cube_derivative(double x) {
    if (x != 0) {
        return -3.0 / (x * x * x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = 1/x^4
double reciprocal_x4_derivative(double x) {
    if (x != 0) {
        return -4.0 / (x * x * x * x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = 1/x^5
double reciprocal_x5_derivative(double x) {
    if (x != 0) {
        return -5.0 / (x * x * x * x * x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = sqrt(x)
double sqrt_func_derivative(double x) {
    if (x > 0) {
        return 0.5 / sqrt(x);
    } else {
        return 0.0;
    }
}

// Derivative of y = sqrt(x)^3
double cube_root_derivative(double x) {
    if (x > 0) {
        return 1.5 * sqrt(x);
    } else {
        return 0.0;
    }
}

// Derivative of y = 1/sqrt(x)
double reciprocal_sqrt_derivative(double x) {
    if (x > 0) {
        return -0.5 / (x * sqrt(x));
    } else {
        return 0.0;
    }
}

// Derivative of y = exp(x)
double exp_func_derivative(double x) {
    return exp(x);
}

// Derivative of y = log(x)
double log_func_derivative(double x) {
    if (x > 0) {
        return 1.0 / x;
    } else {
        return 0.0;
    }
}

// Derivative of y = abs(x)
double abs_func_derivative(double x) {
    if (x > 0) {
        return 1.0;
    } else if (x < 0) {
        return -1.0;
    } else {
        // Technically undefined at x = 0
        return 0.0;
    }
}

// Derivative of y = sin(x)
double sin_func_derivative(double x) {
    return cos(x);
}

// Derivative of y = cos(x)
double cos_func_derivative(double x) {
    return -sin(x);
}

// Derivative of y = tan(x)
double tan_func_derivative(double x) {
    return 1.0 + tan(x) * tan(x);
}

// Derivative of y = tanh(x)
double tanh_func_derivative(double x) {
    return 1.0 - tanh(x) * tanh(x);
}

// Derivative of y = sgn(x)
double sgn_func_derivative(double x) {
    // The derivative of sgn(x) is not defined for all x, 
    // so returning 0 by default, but it is technically a piecewise function.
    return 0.0;
}

// Derivative of y = arcsin(x)
double arcsin_func_derivative(double x) {
    if (x > -1 && x < 1) {
        return 1.0 / sqrt(1.0 - x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = arccos(x)
double arccos_func_derivative(double x) {
    if (x > -1 && x < 1) {
        return -1.0 / sqrt(1.0 - x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = arctan(x)
double arctan_func_derivative(double x) {
    return 1.0 / (1.0 + x * x);
}

// Derivative of y = arctanh(x)
double arctanh_func_derivative(double x) {
    if (x > -1 && x < 1) {
        return 1.0 / (1.0 - x * x);
    } else {
        return 0.0;
    }
}

// Derivative of y = 0
double zero_func_derivative(double x) {
    return 0.0;
}

// Derivative of y = exp(-x^2), often called the Gaussian function
double gauss_func_derivative(double x) {
    return -2.0 * x * exp(-x * x);
}







// 数据生成
void generate_data(double *inputs, double *outputs, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        inputs[i] = ((double)rand() / RAND_MAX) * 6 - 3; // -3到3之间的随机数
        outputs[i] = 3 *exp(0.5 *inputs[i] + 1)+4; // 
    }
}

double calculate_rmse(const double *actual, const double *predicted, int size) {
    double sum_squared_error = 0.0;
    
    for (int i = 0; i < size; i++) {
        double error = actual[i] - predicted[i];
        sum_squared_error += error * error;
    }

    return sqrt(sum_squared_error / size);
}

void standardize_data(double *data, int num_samples) {
    double sum = 0.0;
    double mean = 0.0;
    double variance = 0.0;
    double stddev = 0.0;

    // 计算均值
    for (int i = 0; i < num_samples; i++) {
        sum += data[i];
    }
    mean = sum / num_samples;

    // 计算方差
    sum = 0.0;
    for (int i = 0; i < num_samples; i++) {
        sum += (data[i] - mean) * (data[i] - mean);
    }
    variance = sum / num_samples;

    // 计算标准差
    stddev = sqrt(variance);

    // 标准化数据
    for (int i = 0; i < num_samples; i++) {
        if (stddev != 0) {
            data[i] = (data[i] - mean) / stddev;
        }
    }
}

double normal_random(double sigma) {
    static const double epsilon = 1e-7;
    static const double two_pi = 2.0 * M_PI;

    double u1, u2;
    do {
        u1 = rand() / (double)RAND_MAX;
        u2 = rand() / (double)RAND_MAX;
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);

    return z0 * sigma; // 均值为0，标准差为sigma
}