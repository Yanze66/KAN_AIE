#include <adf.h>
// #include "adf/new_frontend/adf.h"
#include "kernels.h"

using namespace adf;

class simpleGraph : public adf::graph {
private:
//   kernel repeat_node;
  kernel first;     

public:
  input_gmio in;

  output_plio out0;
  output_gmio out;
  simpleGraph()
  {
    //in0 = input_plio::create(plio_32_bits, "data/input_1.txt");
    //in1 = input_plio::create(plio_32_bits, "data/input_1.txt");

    //out0 = output_plio::create(plio_32_bits, "data/output_1.txt");
    //out1 = output_plio::create(plio_32_bits, "data/output_2.txt");
    // repeat_node = kernel::create(repeat);
    first = kernel::create(spline);
    in = input_gmio::create("in",256,2048); //burst_length
    out = output_gmio::create("out",256,2048);

    out0 = output_plio::create(plio_32_bits, "data/output_1.txt");
    //connect(in0.out[0], first.in[0]); 
    //dimensions(first.in[0]) = { 8 };
    //single_buffer(first.in[0]);               //uncomment for single buffer, by default double buffer will be used.
    //connect(first.out[0], out0.in[0]);
    //dimensions(first.out[0]) = { 8 };

    //lut0 = parameter::array (LUT);
    //connect<> (lut0, first);

    
    connect<stream> net0 (in.out[0], first.in[0]);
    connect<stream> net2 (first.out[0], out0.in[0]); //to plio
    connect<stream> net3 (first.out[0], out.in[0]); //to plio
    fifo_depth(net0) = 512;
    //connect<parameter> (vectorInput, first.in[2]);    //connection for RTP_array
    //connect<parameter> (sync(first.inout[0]), vectorOutput);
    
    //connect<parameter> (mlen, first.in[1]);
    //connect<parameter> (new_len, first.in[2]);
    
    // source(repeat_node)="src/kernels/repeat_node.cc";
    source(first) = "src/kernels/neural.cc";
    // runtime<ratio>(repeat_node) = 0.9;
    runtime<ratio>(first) = 0.9;
    }
};

