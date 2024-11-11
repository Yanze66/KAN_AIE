
# Introduction

The folders src and src_active are KAN developed from a piecewise linear function, which differs from the KAN disclosed by MIT (in folder bspli-KAN & kan-bspline vectorized). To be honest, it took me some effort to try two different KANs, and I finally chose the official KAN (MIT version from Liu Ziming) because their version can extend KAN to arbitrary depth and width, which means it has the potential to figure out more complex tasks.

I built the C-based KAN in folder bspline-kan. It's a step-by-step implementation, which develops from the spline function to the whole network. 

In folder kan-bspline-vectorized, I am gonna optimize KAN. So the first step is to vectorize it, which means the core function, ie, Bspline, processes data in a matrix format rather than a singular value.

# ==== Splitline ====Update in 2024.11.11========

# Reference

This repository (src) is developed from http://openkan.org/KANscore.html

Having learned basic knowledge of KAN, I wrote a C library to execute as a benchmark. 


Another folder (src_activate) is developed according https://arxiv.org/pdf/2404.19756

They suggested that by adding an activation, KAN overcomes the COD. 

# HOW TO USE

1. cd KAN_AIE/KAN_c_src/src
   
2. gcc *.c -g -lm -o main 

3. ./main 


# Future work

It‘s based on the piecewise linear function and only two layers. To extend the KAN to arbitrary width and depth, I gonna move to Liu's KAN architecture. 
