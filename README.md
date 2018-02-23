# 1-Bit Tensor Completion

This repository contains MATLAB codes and scripts designed for the 1-bit tensor completion algorithm.
In the proposed approach is explored the recovery of a low-rank tensor from a small number of binary
measurements. Specifically, given a 3-order tensor where only a small number of binary entries are
available, we unfold it into 3 matrices and we apply the quantized matrix completion algorithm to the
all-mode matricizations of the tensor.

## Requirements

### Dataset
The performance of the proposed scheme, is quantified using hyperspectral Earth Observation images
taken from airbornes or satellites which are publicly available in 
http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes
Specifically, we considered the hyperspectral images over Indian Pines for this particular experiment.

### Matrix Completion Algorithm 
We use the matrix completion algorithm, which is available in
http://perception.csl.illinois.edu/matrix-rank/sample_code.html and solves the matrix completion problem
using the Augmented Lagrange Multipliers (ALM) method.

## Contents
demo_1btc.m : The primary script that loads the data, performs the recovery of a low-rank tensor from a
number of binary measurements and provides the results.

**QMC.m**: Perform the quantized matrix completion for each unfolding of the tensor.

**quantization.m**: Quantize to a single bit the entries of the original tensor, using as threshold the
mean value of each pixel.

**bin_boundaries.m**: Compute the upper and lower quantization bin boundaries of each measurement.

**Unfold.m**: Unfold the tensor into a matricization mode.

**Fold.m**: Fold a matricization mode into a tensor.

**f_log.m**: Compute the inverse logit link function.

**f_log_prime.m**: Compute the derivative of the inverse logit link function.

**f_pro.m**: Compute the cdf of a Gaussian distribution at specified values (inverse probit link function).

**f_pro_prime.m**: Compute the derivative of the inverse probit link function.

## References
1. A.  Aidini,  G.  Tsagkatakis,  and P.  Tsakalides,
“1-Bit  Tensor  Completion,” in Proc.  2018 IS&T  International Symposium  on  Electronic  Imaging,
Image  Processing:  Algorithms  and  Systems, Burlingame,  CA,  January  28-February 2, 2018.

2. Lin, Zhouchen, Minming Chen, and Yi Ma.
"The augmented lagrange multiplier method for exact recovery of corrupted low-rank matrices."
arXiv preprint arXiv:1009.5055 (2010).
