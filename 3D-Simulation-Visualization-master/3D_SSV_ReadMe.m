% See Project: Harmonic Analysis of Spherical distributions by 
%%  https://www.researchgate.net/project/Harmonic-Analysis-of-Spherical-distributions
%% Product Overview
% 3D-Directional Statistics, Simulation and Visualization (3D-SSV) is a
% fast, accurate, and convenient way to calculate, simulate, and visualize
% from 3-dimensional spherical distributions.  In addition to the ability
% to simulate from the entire family of Fisher-Bingham distributions, what
% makes this package unique is fourfold.
%
% First, this package makes use of the Heirarchical Equal Area isoLatitude
% Pixelization (HEALPix) of data on the sphere, developed by the Jet
% Propulsion Laboratory
% <http://healpix.jpl.nasa.gov/ HEALPix>
% Additionally, this package also requires particular functions from the
% MEALPix package (Copyright 2010-2011 Lee Samuel Finn) to make use of the
% many desirable characteristics from the HEALPix tesselation.
%
% Additionally, this package can use the spherical harmonic characterization
% of a probability density function to generate random samples and plot both
% the density function as well as the simulated data points.  Included in this
% package, there is one specific order (L = 3) and degree (m = 2) implemented
% for both real valued spherical harmonics and complex valued spherical
% harmonics [see Random_Y3_2Real_square.m and Random_Y3_2Compl_square.m,
% respectively]. Additionally, we developed code that will calculate and plot
% the densities for alternate specifications of degree and order [see
% Density_SphHarm.m and Density_SphHarm_All.m], and we encourage the user to
% explore this area if she or he is so inclined.
%
% Further aspects of the package that distinguish it from others are the plots
% themselves.  When we were looking for appropriate graphical representations
% of the simulated data from these distributions, we could not find anything
% that met our specific needs, so we developed our own.  Contrary to the common
% 3-dimensional heatmap style plots (which are useful, but not what we felt
% would give the best representation of the data in this particular instance),
% 3D-SSV plots the 3-dimensional spherical density plots topographically on
% the surface of the unit sphere, where the values are given by the mathematical
% formulation of the density function for a given set of parameter values (all of
% the 'Density_XXXXX.m' functions will plot the respective density). In addition,
% we developed a method to render histograms on the surface of the unit sphere for
% data points [see Plot_Hist2DSphere.m], regardless of whether that data is
% simulated from the aforementioned densities or collected via alternate
% methods.  This last contribution may be particularly useful to users who wish
% to visualize 3-dimensional directional data in a way that is intuitive, and
% up to this point, largely underutilized (with respect to Matlab implementation).
%
%
% Below, we give some examples that the user may copy, paste, (uncomment) and run
% to get an idea of some of the functionality of this package.
%
%
%
% scrsz = get(groot,'ScreenSize');
% set(0,'defaultfigureposition',[scrsz(3)-650 scrsz(4)-650 600 550]')
%
% clear all
% close all
% n=2^13;
% nSide=2^3;
% nPix= nSide2nPix(nSide);

% kappa = 4.2;
% bet = 4.5 ;  % beta
% gamm = -3.5; % gamma
% Psi=  pi/2;  °% for rotation to the new frame of reference
%
% Mu=[0 0 1];  % for rotating the North pole to Mu
%
% ensure unit length, also done internally within the appropriate functions
% Mu=Mu/norm(Mu);
%
% resolution=100;  % graphical parameter
%
%
%
%% Uniform Distribution on S2 (unit sphere)
%
% generate n random points from the uniform distribution on S2
% Y = Random_Uni_Inv(n);
% hY = Hist2DSphere(Y,nPix);
%
% Uniform Density Plot
% figure;
% Plot_DataRandomS2(Y);
% title(['Uniform; Mean =[', num2str(mean(Y)),']'])% /norm(mean(Y)
%
% Uniform Histogram Plot
% figure;
% Plot_Hist2DSphere(hY);
% title(['Uniform; N=', num2str(n),'; nPix=', num2str(nPix)])
%
%
%
%% von Mises Fisher
%
%
%% von Mises Fisher Density Plot
% figure;
% gx1 = Density_vMF(kappa,Mu,resolution);
% Y = Random_vMF_Inv(kappa,Mu,n);
% hY=Hist2DSphere(Y,nPix);
%
%% von Mises Fisher simulated data points plot
% figure;
% hp=Plot_DataRandomS2(Y);
% title(['vMF; Mean =[', num2str(mean(Y)),']'])% /norm(mean(Y)
%
%% von Mises Fisher Histogram Plot
% figure;
% hh=Plot_Hist2DSphere(hY);
% title(['von Mises Fisher; \kappa =', num2str(kappa),'; Max=', num2str(max(hY))])
%
%
%% Fisher-Bingham 4: Density + Random Sample + Histogram
% kappa=0;
% figure;
% gx4 = Density_FB4(kappa,gamm,Mu,resolution);
                                                                                                                                                                                                                                        
% Y = Random_FB4(kappa,gamm,Mu,n);
% figure;
% Plot_DataRandomS2(Y);
% title(['FB4; Mean =[', num2str(mean(Y)),']'])% /norm(mean(Y)

% hY=Hist2DSphere(Y,nPix);
% figure;
% Plot_Hist2DSphere(hY);
% title(['FB4; \kappa =', num2str(kappa), '; \gamma =',num2str(gamm)]);
%
%% This work of Gy.Terdik is supported by the 
% EFOP-3.6.2-16-2017-00015 project. The project has been
% supported  by the European Union, co-financed by the European 
% Social Fund.
%
%
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik
% gyorgy.terdik@gmail.com


