function Y = Random_FB4_Beta_GyT(bet,varargin)
% RANDOM_BETA_FB4_GYT Generates a random sample of size N from the Fisher-Bingham_5 distribution
%
% Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% bet = 4.5 ;
% Psi= pi/2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n)
% Y = Random_FB4_Beta_GyT(bet,Mu,Psi)
% Y = Random_FB4_Beta_GyT(bet,Mu)
% Y = Random_FB4_Beta_GyT(bet)
%%
%
% Required Inputs
% bet         distribution parameter (beta >= 0)
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
noneg = {'nonnegative'};
class = {'numeric'};
pos_int = {'integer','positive'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'bet',@(x)validateattributes(x,dbl,noneg));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(bet,varargin{:});
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
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2017 Gy.Terdik, B.Wainwright
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
%% 

if  bet==0
    Y = Random_Uni_Norm(n);
    return
end 

%% density
% betT = acosh(besseli(0,bet));
C_minus = integral(@(u)(exp(-bet*u.^2)),-1,1);
C_plus = integral(@(u)(exp(bet*u.^2)),-1,1);
C = exp(bet)*C_minus+exp(-bet)*C_plus;
p1 = exp(bet)*C_minus/(exp(bet)*C_minus+exp(-bet)*C_plus);
X=[];
while length(X)<n
    U = rand(n,1);
    V1_minus = Watson_LW(-bet,n);
    V1_plus = Watson_LW(bet,n);
    V = (U<p1).*V1_minus+(1-(U<p1)).*V1_plus;
    X2 = V(U <= besseli(0,bet*(1-V.^2))./cosh(bet*(1-V.^2)));
    X = [X;X2];
end
X = X(1:n);

%% phi calculating the x1 and x2 components;
Phi = vMF_Circ_Phi(bet*(1-X.^2));  % vMF
% U = randi(4,n,1);
% Phi = acos(PhivMF)/2; %+((U<=1/2)-(U>1/2))*pi
% Phi = Phi.*((U==1)+(U==3))+ ((U==2)+(U==3))*pi-Phi.*((U==2)+(U==4))+(U==4)*2*pi; % doing 2 pi
cPhi = cos(Phi+Psi);
sPhi = sin(Phi+Psi);

%%
sX = sqrt(1-X.^2); %
Y = [cPhi.*sX,sPhi.*sX,X];

%% rotation (orient the North Pole with mean vector (mu))
Np=[0,0,1]; % z-axis (North Pole)
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end