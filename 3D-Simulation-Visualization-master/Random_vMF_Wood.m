function Y = Random_vMF_Wood(kappa,varargin)
% RANDOM_VMF_WOOD Generates a random sample of size N from the 
%                 von Mises-Fisher distribution
%
% Y = Random_vMF_Wood(kappa,Mu,n,d);
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_vMF_Wood(kappa,Mu,n,d);
% Y = Random_vMF_Wood(kappa,Mu,n);
% Y = Random_vMF_Wood(kappa,Mu);
% Y = Random_vMF_Wood(kappa)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter, (for rotating the North pole to Mu)
% n             sample size ( where n is a positive integer)
% d             dimension of the space, for a 2D sphere d=3

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
class = {'numeric'};
pos_int = {'integer','positive'};
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));
addOptional(p,'d',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;
d = p.Results.d;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^15;
end
if isempty(d)
    d = 3;
end

%% 
% d=3;   dim of the space / m=3 corresponds to the unit sphere

if kappa==0;
    Y = RandUni_Norm_N(n); % kappa=0 corresponds to uniform dist. in S2
else
    W1 = vMF_Wood(kappa,n,d);
end

%%
V = randn((d-1),n);
V = V./kron(ones(d-1,1),sqrt(sum(V.^2,1)));
Y = [kron(ones(1,d-1),(1-W1.^2).^(1/2)).*V' W1];
%% 
if d  ~= 3
    return
end

%% rotation (orient mean vector (mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
return