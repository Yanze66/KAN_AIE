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
