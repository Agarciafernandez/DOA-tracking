This repository contains the Matlab implementations of the (sigma-point) iterated posterior linearisation filter (IPLF) for single-target tracking with direction-of-arrival (DOA) measurements. The filter consisders von Mises-Fisher DOA measurements [1] and Kent DOA measurements [2]. The filters obtain a Gaussian approximation of the posterior based on iterated statistical linear regressions in the update step. In the code, the DOA measurement data for both filters can be chosen to be von Mises-Fixher/Kent or Gaussian distributed. 


This repository makes use of the 3D-Directional-SSV package to sample the Kent distribution: 
https://github.com/TerdikGyorgy/3D-Simulation-Visualization.


The explanation of the IPLF can be found in [3]-[4]. The work [2] was developed as part of the DSTL Grant no. 1000143726.


- IPLF_3D_VMF_filter.m runs the IPLF with von Mises Fisher DOA measurements.

- IPLF_3D_Kent_filter.m runs the IPLF with Kent DOA measurements.

References

[1] Á. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises–Fisher Measurements," in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, 1 June1, 2019

[2] A. F. García-Fernández, S. Maskell, P. Horridge, J. Ralph, “Gaussian tracking with Kent-distributed direction-of-arrival measurements,” IEEE Transactions on Vehicular Technology, 2021.

[3] Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015

[4] F. Tronarp, Á. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing in Nonlinear and Non-Gaussian Systems Using Conditional Moments," in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018
