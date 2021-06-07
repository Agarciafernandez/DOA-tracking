function Y = Random_Watson_Fish(gamm,varargin)
% RANDOM_WATSOM_FISH Generates a random sample of size n from the 
%                    Watson distribution
%
% Y = Random_Watsom_Fish(gamm,Mu,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% gamm = -3.5; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_Watsom_Fish(gamm,Mu,n)
% Y = Random_Watsom_Fish(gamm,Mu)
% Y = Random_Watsom_Fish(gamm)
%%
%
% Required Inputs
% gamm        distribution parameter (gamma \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter, (for rotating the North pole to Mu)
% n           sample size ( where n is a positive integer)
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
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'gamm',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(gamm,varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^12;
end


%% when gamma==0 Watson_Fish is uniform distribution
if gamm==0;
    Y = Random_Uni_Norm(n);
    return
else 
    X = Watson_Fish(gamm,n);  
end

%% psi iid uniform
psi = 2*pi*rand(length(X),1);
sX=sqrt(1-X.^2);
Y=[cos(psi).*sX,sin(psi).*sX,X];

%% rotation (orient mean vector (Mu) with the North Pole)
Mu = Mu/norm(Mu); % should be unit vector
Np=[0,0,1]; % z-axis (North Pole)
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end