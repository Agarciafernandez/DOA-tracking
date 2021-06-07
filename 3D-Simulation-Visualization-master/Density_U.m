function gx = Density_U(varargin)
% DENSITY_U Calculates and plots the U-density (not to be confused with the
% uniform density)
%
% gx = Density_U(Mu,Psi, resolution);
%%
% Examples of correct usage:
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
% resolution=100;
%
% Examples of correct function construction:
% gx = Density_U(Mu, resolution)
% gx = Density_U(Mu) 
% gx = Density_U
%%
%
% Optional Inputs
% Mu          distribution parameter, (Mu:mean vector with finite elements)         
% resolution  plot resolution parameter
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik
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
vec3D = {'size',[1,3],'finite'};
class = {'numeric'};
pos_int = {'integer','positive'};
longPsi = {'>=',0,'<=',2*pi,'finite'};

% Check characteristics for the required and optional parameters
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,vec3D,'Mu'));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,longPsi,'Psi'));
addOptional(p,'resolution',[],@(x)validateattributes(x,class,pos_int));

p.parse(varargin{:});
Mu = p.Results.Mu;
Psi = p.Results.Psi;
resolution = p.Results.resolution;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(Psi)
    Psi = 0;
end
if isempty(resolution)
    resolution = 100;
end

%% 

delta= pi/resolution; % step size
theta = 0:delta/2:pi/2; % colatitude-delta
phi = 0:delta: pi/4; % longitude-delta

%% plot
[phi1,theta1] = meshgrid(phi,theta);

Phi = [phi, (phi+pi/4), (phi+pi/2), (phi+3*pi/4),phi+pi, (phi+5*pi/4), (phi+3*pi/2), (phi+7*pi/4)]; %  0:delta:2*pi -delta; %  -delta
Phi = Phi-Psi;
Theta= [theta,theta+pi/2 ] ;%  0:delta/2:pi ;
[Phi1,Theta1] = meshgrid(Phi,Theta);
%% 

x2 = sin(Theta1).*cos(Phi1);
y2 = sin(Theta1).*sin(Phi1);
z2 = cos(Theta1);
%% 
tg_th=tan(theta);
cos_phi1=1./cos(phi);
[cos_1,tg_1]=meshgrid(cos_phi1,tg_th);
ind_cos=[tg_1<cos_1];
sinTheta=sin(theta1);
sinTheta(sinTheta==0)=1;
gx1=cos(phi1)./cos(theta1).^2./cos(phi1).^2/16.*ind_cos;
gx2=1./sinTheta.^2./cos(phi1).^2.*cos(phi1)/16.*(1-ind_cos);
gx14=gx1+gx2 ;
gx12=[gx14 fliplr(gx14)];
Gx1=[gx12 gx12 gx12 gx12];
gx=[Gx1; flipud(Gx1)];

%% 


r = gx;
x = (1+r).*x2;
y = (1+r).*y2;
z = (1+r).*z2;
h = surf(x,y,z,r+1); 

%% rotation (orient mean vector (Mu) with the North Pole)
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
set(h,'LineStyle','none')
colormap(jet(1024));
colorbar
camlight left
camlight right
lighting phong
A=max(1.5,max(1+gx(:)));
axis([-A A -A A -A A])
grid on

% title('Uniform Distribution')
xlabel('x-axis')
ylabel('y-axis')