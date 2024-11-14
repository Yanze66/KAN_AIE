#ifndef spline_h
#define spline_h

// #include "adf/x86sim/streamStructs.h"
#include <adf.h>

// void spline_vec(input_stream<float> * __restrict bufin, output_stream<float>* __restrict bufout);
void spline1(input_stream<float> *  __restrict datain, output_stream<float>* dataout, output_stream<float>* resultout);
void spline2(input_stream<float> *  __restrict datain, output_stream<float>* dataout, input_stream<float> *  __restrict resultin,output_stream<float>* resultout);
void spline3(input_stream<float> *  __restrict datain,  input_stream<float> *  __restrict resultin,output_stream<float>* resultout);
// void split(input_stream<float> *  __restrict datain,  output_stream<float>* dataout1,output_stream<float>* dataout2);
// void merge(input_stream<float> *  __restrict datain1, input_stream<float> *  __restrict datain2, output_stream<float>* dataout);

// void repeat(input_stream<float> * __restrict bufin, output_stream<float>* __restrict bufout);




#endif
