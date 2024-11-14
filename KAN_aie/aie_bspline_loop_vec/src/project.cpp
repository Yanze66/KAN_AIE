#include "project.h"
// #include <cmath>
#include <fstream>
#include <cstdlib>
#include <string.h>
// #include <inttypes.h>
// #include <sys/time.h>
#include "math.h"
simpleGraph mygraph;


#if defined(__AIESIM__) || defined(__X86SIM__)
// static unsigned long long cpucycles(void)
// {
//     unsigned long long result;
//       __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
//               : "=a" (result) ::  "%rdx");
//         return result;
// }

// 数据生成
void generate_data(float *inputs, float *outputs, int num_samples) {
    for (int i = 0; i < num_samples; i++) {
        inputs[i] = ((float)rand() / RAND_MAX) * 6 - 3; // -3到3之间的随机数
        outputs[i] = 1.5 *exp(0.5 *inputs[i] + 1)+1; // 
    }
}

int main(int argc, char **argv) 
{
    // uint64_t cyc0; 
    // uint64_t cyc1;
    // int num_samples=1000;
    // float inputs[num_samples];
    // float outputs[num_samples];
    // generate_data(inputs, outputs, num_samples);

    float inputs[30]={-0.5,0,0.5};
    
    float  *inputbuf  = (float*) GMIO::malloc(30*sizeof(float)); //GMIO input in this way
    float *outputbuf = (float*) GMIO::malloc(24*sizeof(float));

    

    printf("\nThe input to nerural is :\n");
    	for(int i=0;i<30;i++)
	{
        inputbuf[i]=inputs[i];
        // inputbuf[2*i+1]=outputs[i];
		printf("%f",inputbuf[i]);	
	}

    // cyc0 = cpucycles();
    mygraph.init();
    
    
    mygraph.run(1);
    for(int i=0;i<1;i++){
        mygraph.in.gm2aie_nb(&inputbuf[0],4*sizeof(float));
        //mygraph.update(mygraph.mlen, mlen);
  	    //mygraph.wait(); 
        mygraph.out.aie2gm_nb(&outputbuf[0],24*sizeof(float));
        mygraph.out.wait(); //important !! wait for data move to DDR
    }
    mygraph.end();

    // cyc1 =  cpucycles();
    // printf("\ncycles of sha256 : %" PRIu64 "\n\n", cyc1-cyc0);
    
    printf("\nThe output is :\n");  //  a converse happens,ie, ABCD-> DCBA, so we print in conversed order
    for(int i=0;i<24;i++){	
        printf("%02f",outputbuf[i]);
                   	
	}
    printf("\n");

    GMIO::free(inputbuf);
    GMIO::free(outputbuf);
  return 0;
} 
#endif