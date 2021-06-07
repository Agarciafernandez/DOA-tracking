This repository contains the Matlab implementations of the (sigma-point) iterated posterior linearisation filter (IPLF) with direction-of-arrival (DOA) measurements with von Mises-Fisher distribution [1] and Kent distribution [2]. 

The DOA measurement data for both filters can be chosen to be Kent or Gaussian distributed. 


This repository makes use of the 3D-Directional-SSV package to sample the Kent distribution: 
https://github.com/TerdikGyorgy/3D-Simulation-Visualization.


The explanation of the IPLF can be found in [3]-[4]. The work [2] was developed as part of the DSTL Grant no. 1000143726.

References

[1] �. F. Garc�a-Fern�ndez, F. Tronarp and S. S�rkk�, "Gaussian Target Tracking With Direction-of-Arrival von Mises�Fisher Measurements," in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, 1 June1, 2019

[2] A. F. Garc�a-Fern�ndez, S. Maskell, P. Horridge, J. Ralph, �Gaussian tracking with Kent-distributed direction-of-arrival measurements,� IEEE Transactions on Vehicular Technology, 2021.

[3]�. F. Garc�a-Fern�ndez, L. Svensson, M. R. Morelande and S. S�rkk�, "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015

[4]F. Tronarp, �. F. Garc�a-Fern�ndez and S. S�rkk�, "Iterative Filtering and Smoothing in Nonlinear and Non-Gaussian Systems Using Conditional Moments," in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018
