function gx = Density_vMF(kappa,varargin)
% DENSITY_VMF Calculates and plots von Mises-Fisher density
%
% gx = Density_vMF(kappa,Mu,resolution);
%%
% Examples of correct usage:
%
% kappa = 4.2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
% resolution=100;
%
% Examples of correct function construction: 
% gx = Density_vMF(kappa,Mu,resolution) 
% gx = Density_vMF(kappa,Mu)
% gx = Density_vMF(kappa)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter, (Mu:mean vector with finite elements)          
% resolution  plot resolution parameter
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik, B.Wainwright
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
dbl = {'double'};
class = {'numeric'};
pos_int = {'integer','positive'};
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'resolution',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,varargin{:});
Mu = p.Results.Mu;
resolution = p.Results.resolution;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(resolution)
    resolution = 100;
end

%% density
cFB=(exp(kappa)-exp(-kappa))/kappa; % normalization constant 
delta = pi/resolution; % step size
theta = 0:delta:pi; % colatitude
phi = 0:2*delta:2*pi; % longitude

%% plot
[phi,theta] = meshgrid(phi,theta);
gx = exp(kappa*cos(theta))/cFB/2/pi; % exp(kappa)/cFB/2/pi
x2 = sin(theta).*cos(phi);
y2 = sin(theta).*sin(phi);
z2 = cos(theta);
                 
r = gx;
x = (1+r).* x2;
y = (1+r).*y2;
z = (1+r).*z2;
h = surf(x,y,z,r+1);

%% rotation (orient mean vector (mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if Mu(3)~= 1
    Ux=cross(Np,Mu);%axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=rad2deg(acos(Mu(3)));
    rotate(h,Ux,thetaX );
end

%% display options
view(3)
set(h,'LineStyle','none');
colormap(jet(1024));
colorbar
camlight left
camlight right
lighting phong
A=max(1.5,max(1+gx(:)));
axis([-A A -A A -A A])
grid on
                 
title(['von Mises Fisher; \kappa =', num2str(kappa),' Max=',num2str(max(gx(:)))])
