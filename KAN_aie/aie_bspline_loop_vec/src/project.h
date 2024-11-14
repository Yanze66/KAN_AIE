#include <adf.h>
// #include "adf/new_frontend/adf.h"
#include "adf/new_frontend/adf.h"
#include "kernels.h"

using namespace adf;

class simpleGraph : public adf::graph {
private:
//   kernel repeat_node;
  kernel first;     
  kernel second;
  kernel third; 

//   kernel split1;
//   kernel split2;

//   kernel merge1;
//   kernel merge2; 
public:
  input_gmio in;

  output_plio out0;
  output_gmio out;
  simpleGraph()
  {
   
    first = kernel::create(spline1);
    second = kernel::create(spline2);
    third = kernel::create(spline3);

    // split1 = kernel::create(split);
    // split2 = kernel::create(split);

    // merge1 = kernel::create(merge);
    // merge2 = kernel::create(merge);

    in = input_gmio::create("in",128,256); //burst_length
    out = output_gmio::create("out",128,256);

    out0 = output_plio::create(plio_32_bits, "data/output_1.txt");
   

    
    // connect<stream> net0 (in.out[0], split1.in[0]);
    // connect<stream> net1 (split1.out[0], split2.in[0]); //to plio
    // connect<stream> net2 (split1.out[1], third.in[0]);



    // connect<stream> net3 (split2.out[0], first.in[0]); //to plio
    // connect<stream> net4 (split2.out[1], second.in[0]);

    // connect<stream> net5 (first.out[0], merge1.in[0]); //to plio
    // connect<stream> net6 (second.out[0], merge1.in[1]); //to plio

    // connect<stream> net7 (merge1.out[0], merge2.in[0]); //to plio
    // connect<stream> net8 (third.out[0], merge2.in[1]); //to plio

    // connect<stream> net9 (merge2.out[0], out.in[0]);
    // connect<stream> net10 (merge2.out[0], out0.in[0]);

    connect<stream> net0 (in.out[0], first.in[0]);

    connect<stream> net2 (first.out[0], second.in[0]);
    connect<stream> net3 (second.out[0], third.in[0]);

    connect<stream> net1 (first.out[1], second.in[1]);
    connect<stream> net4 (second.out[1], third.in[1]);
    connect<stream> net5 (third.out[0], out0.in[0]);
    connect<stream> net6 (third.out[0], out.in[0]);

    fifo_depth(net0) = 128;
    // fifo_depth(net1) = 512;
    fifo_depth(net2) = 512;
    fifo_depth(net3) = 512;
    fifo_depth(net4) = 512;
    fifo_depth(net5) = 512;
    fifo_depth(net6) = 512;
    // fifo_depth(net7) = 512;
    // fifo_depth(net8) = 512;
    // fifo_depth(net9) = 512;
    // fifo_depth(net10) = 512;

    
    source(first) = "src/kernels/neural.cc";
    source(second) = "src/kernels/neural.cc";
    source(third) = "src/kernels/neural.cc";

    // source(split1)    = "src/kernels/merge_split.cc";
    // source(split2)    = "src/kernels/merge_split.cc";

    // source(merge1)    = "src/kernels/merge_split.cc";
    // source(merge2)    = "src/kernels/merge_split.cc";
    
    // runtime<ratio>(repeat_node) = 0.9;


    runtime<ratio>(first) = 0.1;
    runtime<ratio>(second) = 0.1;
    runtime<ratio>(third) = 0.1;

    // runtime<ratio>(split1) = 0.1;
    // runtime<ratio>(split2) = 0.1;
    // runtime<ratio>(merge1) = 0.1;
    // runtime<ratio>(merge2) = 0.1;
    adf::location<kernel>(first)=adf::tile(26,0); 
    adf::location<kernel>(second)=adf::tile(27,0); 
    adf::location<kernel>(third)=adf::tile(28,0); 
    }
};

