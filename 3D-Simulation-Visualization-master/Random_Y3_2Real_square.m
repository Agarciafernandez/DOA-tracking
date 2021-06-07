function Y = Random_Y3_2Real_square(varargin)
% RANDOM_Y3_2REAL_SQUARE generates a random sample of size n from the
%                        distribution characterized by the real spherical
%                        harmonic of degree l=3 and order m=2
%
% Y = Random_Y3_2Real_square(Mu,Psi,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% Psi= pi/2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_Y3_2Real_square(Mu,Psi,n)
% Y = Random_Y3_2Real_square(Mu,Psi)
% Y = Random_Y3_2Real_square(Mu)
% Y = Random_Y3_2Real_square
%%
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
% Psi         distribution parameter, (-infinity < Psi < infinity,
%                               for rotation to the new frame of reference)
% n           sample size (where n is a positive integer)
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

%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
dbl = {'double'};
class = {'numeric'};
pos_int = {'integer','positive'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(varargin{:});
Mu = p.Results.Mu;
Psi = p.Results.Psi;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(Psi)
    Psi = 0;
end
if isempty(n)
    n = 2^12;
end

%% spherical harmonics parameters: fix degree and order
ld=3;  % degree of spherical harmonic
mo=2;  % order of spherical harmonic

%% W, as simulated in 4.1 (Density 2g_3^2) via acceptance-rejection sampling

%% envelope parameters
% X~e(x) = Beta(alpha=3.08,beta=2.5249) - Envelope Distribution
% U~Unif(0,1)
% Reject X if U > f(X)/(M*e(X)) - Rejection Criteria

M = 1.074; % there should be an avg of M trials per accepted X
alpha = 3.08; 
beta = 2.5249;
n = floor(M*n); % total number of trials necessary

%%  Step 1: Simulate theta 

X2 = betarnd(alpha,beta,n,1); % random n-by-1 vector from e(x)

% Squared, fully normalized associated legendre function of degree ld and 
% orders mo=0,1,...,ld
glm = (legendre(ld,X2,'norm')).^2; 
X1 = 2*glm(mo+1,:);

% envelope density function
g_X2 = betapdf(X2,alpha,beta);

% random n-by-1 vector from uniform(0,1)
U1 = rand(n,1); 
% accepted observations
X_accept = X2(U1<=X1'./g_X2); 

% Positive or Negative w/ Prob(1/2)
U2 = rand(length(X_accept),1); %random vector from uniform(0,1)
W=(ones(length(X_accept),1)-2*(U2<0.5)).*X_accept;

%  theta = arccos(W)
theta = acos(W);

%% Step 2: Simulate phi 

alpha2 = 3/2;
beta2 = 1/2;

LT = length(theta);

X3 = betarnd(alpha2,beta2,LT,1);
U3 = rand(LT,1);
Y1  = zeros(LT,1);
K = 2; % corresponds to cos(ck*x)^2
Nk = K*4;
for k=1:Nk
    s1(k)=sum(and(U3>=(k-1)/Nk, U3< k/Nk)); %count total
    Y1(and(U3>=(k-1)/Nk, U3< k/Nk))= ...
        floor(k/2)*pi+(-1)^(k+1)*acos(X3(and(U3>=(k-1)/Nk,U3<k/Nk)).^(1/2));
end

phi = Y1/K+Psi;

theta_phi = [theta(1:length(phi)) phi];
theta_phi = theta_phi';
theta_phi = mat2cell(theta_phi,2,ones(length(phi),1));
              

%%  Step 3: calculate the point x(theta,phi) from our (theta,phi)
xyz = ang2vec(theta_phi);
A = cell2mat(xyz);
Y = A';

%% rotation (orient the North Pole with mean vector (mu))
Mu = Mu/norm(Mu); % should be unit vector
Np=[0,0,1]; % z-axis (North Pole)
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
