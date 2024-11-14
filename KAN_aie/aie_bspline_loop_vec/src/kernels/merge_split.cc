// #include <stdio.h>
// #include <stdlib.h>
// #include <cstdint>
// #include <string.h>
// #include <cstdio>
// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <adf.h>
// #include <aie_api/aie.hpp>
// #include "aie_api/aie_types.hpp"
#include "../kernels.h"
#include <cstdio>
// #include "adf/x86sim/streamStructs.h"
// #include "adf/x86sim/streamApi.h"

void split(input_stream<float> *  __restrict datain,  output_stream<float>* dataout1,output_stream<float>* dataout2)
{ 
       
    char pipe =3;
    float buf[3];
    for(int p=0;p<pipe;p++){
        buf[p] = readincr(datain);       
    }

  
    // printf("\ninput %d\t",i);
    for(int j=0;j<pipe;j++){           
        writeincr(dataout1, buf[j]);
        writeincr(dataout2, buf[j]);
    }
     
}

void merge(input_stream<float> *  __restrict datain1, input_stream<float> *  __restrict datain2, output_stream<float>* dataout)
{ 
    printf("begin merge")  ; 
    char pipe =7;
    float buf[7*2];
    for(int q=0;q<pipe;q++){
        buf[q] = readincr(datain1); 
    }

    for(int p=0;p<pipe;p++){       
        buf[p+pipe] = readincr(datain2); 
    }  
    // printf("\ninput %d\t",i);
    for(int j=0;j<2*pipe;j++){           
        writeincr(dataout, buf[j]);
        
    }
     
}