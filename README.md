# Simultaneous sparse recovery and blind demodulation
Code to reproduce the figures in the IEEE Transactions on Sigal Processing (TSP) paper "[Simultaneous Sparse Recovery and Blind Demodulation](https://ieeexplore.ieee.org/document/8805114)"

# Abstract
The task of finding a sparse signal decomposition in an overcomplete dictionary is made more complicated when the signal undergoes an unknown modulation (or convolution in the complementary Fourier domain). Such simultaneous sparse recovery and blind demodulation problems appear in many applications including medical imaging, super resolution, self-calibration, etc. In this paper, we consider a more general sparse recovery and blind demodulation problem in which each atom comprising the signal undergoes a distinct modulation process. Under the assumption that the modulating waveforms live in a known common subspace, we employ the lifting technique and recast this problem as the recovery of a column-wise sparse matrix from structured linear measurements. In this framework, we accomplish sparse recovery and blind demodulation simultaneously by minimizing the induced atomic norm, which in this problem corresponds to the block L1 norm minimization. For perfect recovery in the noiseless case, we derive near optimal sample complexity bounds for Gaussian and random Fourier overcomplete dictionaries. We also provide bounds on recovering the column-wise sparse matrix in the noisy case. Numerical simulations illustrate and support our theoretical results.

# Tested on 
Matlab R2017b with [CVX toolbox](http://cvxr.com/cvx/)

# Citation
If you use our method and/or codes, please cite our paper

```
@article{xie2019simultaneous,
  title={Simultaneous sparse recovery and blind demodulation},
  author={Xie, Youye and Wakin, Michael B and Tang, Gongguo},
  journal={IEEE Transactions on Signal Processing},
  volume={67},
  number={19},
  pages={5184--5199},
  year={2019},
  publisher={IEEE}
}
```
