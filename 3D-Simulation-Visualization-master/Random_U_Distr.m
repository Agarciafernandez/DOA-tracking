function Y = Random_U_Distr(varargin)
% RANDOM_U_DISTR Generates a random sample of size n from the 
%                U-distribution (not to be confused with the
%                uniform density)
%
% Y = Random_U_Distr(Mu,Psi,n);
%%
% Examples of correct usage:
%
% n=2^12; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_U_Distr(Mu,n)
% Y = Random_U_Distr(Mu)
% Y = Random_U_Distr
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

%% density 
R1 = rand(3,n)-1/2;
R1U = sqrt(sum(R1.^2,1));
R1U = kron(ones(3,1),R1U);
Y = R1'./R1U';

%% rotation (orient the North Pole with mean vector (mu))
Mu=Mu/norm(Mu); % should be unit vector
Np=[0,0,1]; % z-axis (North Pole)
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
return