#ifndef spline_h
#define spline_h

#include <adf.h>

// len is the len of message, new_len should be pckaged message
void spline(input_stream<float> * __restrict bufin, output_stream<float>* __restrict bufout);
// void repeat(input_stream<float> * __restrict bufin, output_stream<float>* __restrict bufout);




#endif
