function Y = Random_FB4(kappa,gamm,varargin)
% RANDOM_FB4 Generates a random sample of size n from the Fisher-Bingham_4 
%            distribution
%
% Y = Random_FB4(kappa,gamm,Mu,n);
%
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
% gamm = -3.5; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB4(kappa,gamm,Mu,n);
% Y = Random_FB4(kappa,gamm,Mu);
% Y = Random_FB4(kappa,gamm);
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% gamm        distribution parameter (gamma \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
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
addRequired(p,'gamm',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,gamm,varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^12;
end
%% 

%%
sign_kappa=sign(kappa)+(kappa==0);
kappa=abs(kappa); % kappa should be non-negative

%%
if kappa==0 && gamm==0;
    U = rand(n,1);
    X = 2*U-1; % cos(theta)
   
elseif kappa==0 && gamm~=0;
    X= Watson_LW(gamm,n); % Watson !!!!!!!! which one??
    
elseif gamm==0  && kappa ~=0;
    X=vMF_Wood(kappa,n,3);% when gamma==0 FB4 is the von Mises-Fisher distribution
    
else

   X = FB4(kappa,gamm,n); % marginal FB4 distribution ceil(n/a4)

end
%% psi iid uniform
psi = 2*pi*rand(length(X),1);
Xs = sqrt(1-X.^2); 
Y = sign_kappa*[cos(psi).*Xs,sin(psi).*Xs,X];

%% rotation (orient the North Pole with mean vector (Mu))
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end