# 3D-Directional-SSV 
3D-Directional Statistics, Simulation and Visualization (3D-Directional-SSV) is a fast, accurate, and convenient way to calculate, simulate, and visualize from 3-dimensional spherical distributions.  In addition to the ability to simulate from the entire family of Fisher-Bingham distributions, what makes this package unique is fourfold.

First, this package makes use of the Hierarchical Equal Area isoLatitude Pixelization (HEALPix) of data on the sphere, developed by the Jet Propulsion Laboratory\
<http://healpix.jpl.nasa.gov/HEALPix>\
Additionally, this package also requires particular functions from the MEALPix package (Copyright 2010-2011 Lee Samuel Finn) to make use of the many desirable characteristics from the HEALPix tessellation.

Additionally, this package can use the spherical harmonic characterization of a probability density function to generate random samples and plot both the density function as well as the simulated data points.  Included in this package, there is one specific order (L = 3) and degree (m = 2) implemented for both real valued spherical harmonics and complex valued spherical harmonics [see Random_Y3_2Real_square.m and Random_Y3_2Compl_square.m,respectively]. Additionally, we developed code that will calculate and plot
the densities for alternate specifications of degree and order [see Density_SphHarm.m and Density_SphHarm_All.m], and we encourage the user to explore this area if she or he is so inclined.

Further aspects of the package that distinguish it from others are the plots themselves.  When we were looking for appropriate graphical representations of the simulated data from these distributions, we could not find anything that met our specific needs, so we developed our own.  Contrary to the common 3-dimensional heatmap style plots (which are useful, but not what we felt would give the best representation of the data in this particular instance), 3D-SSV plots the 3-dimensional spherical density plots topographically on the surface of the unit sphere, where the values are given by the mathematical formulation of the density function for a given set of parameter values (all of the 'Density_XXXXX.m' functions will plot the respective density). In addition,
we developed a method to render histograms on the surface of the unit sphere for data points [see Plot_Hist2DSphere.m], regardless of whether that data is simulated from the aforementioned densities or collected via alternate methods.  This last contribution may be particularly useful to users who wish to visualize 3-dimensional directional data in a way that is intuitive, and up to this point, largely underutilized (with respect to Matlab implementation).
Copyright 2018 Gy.Terdik
gyorgy.terdik@gmail.com
 See Project: Harmonic Analysis of Spherical distributions by 
 https://www.researchgate.net/project/Harmonic-Analysis-of-Spherical-distributions



