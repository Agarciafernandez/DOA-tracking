function Y = Random_FB5_GyT(kappa,bet,varargin)
% RANDOM_FB5_GYT Generates a random sample of size n from the Fisher-Bingham_5 
%                distribution
%
% Y = Random_FB5_GyT(kappa,bet,Mu,Psi,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
% bet = 4.5 ; 
% Psi= pi/2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB5_GyT(kappa,bet,Mu,Psi,n)
% Y = Random_FB5_GyT(kappa,bet,Mu,Psi)
% Y = Random_FB5_GyT(kappa,bet,Mu)
% Y = Random_FB5_GyT(kappa,bet)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% bet         distribution parameter (beta >= 0)
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
% Psi         distribution parameter, (-infinity < Psi < infinity, 
% %                             for rotation to the new frame of reference)
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
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addRequired(p,'bet',@(x)validateattributes(x,dbl,noneg));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,bet,varargin{:});
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

%% 
if kappa==0 && bet==0
    Y=Random_Uni_Norm(n);
    return
end 
if  bet==0
    Y = Random_vMF_Wood(kappa,Mu,n);
    return
end 

if kappa==0 
   Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n);
    return
end 

%% Density
sign_kappa = sign(kappa); % +(kappa==0)
kappa = abs(kappa); % kappa should be non-negative
X1=[]; X2=[];

% calculating the x3 component of the (x1,x2,x3) coordinate and 
% mixing proportion (p2)
if 2*bet <= kappa % unimodal case
    a1 = kappa-2*bet;
    a2 = kappa+2*bet;
    sig1 = sqrt(1/(4*a1+ sqrt(8*bet)));
    sig2 = sqrt(1/4/kappa);
    gDens = @(A,c0,u) (exp(c0-A*u.^2));%/const_e(lam0)
    kDens = @(a0,lam0,u) (exp(-(2*a0*u.^2+lam0*u.^4)));%/const_B(a0,lam0)
    while(length(X2) < n)
        %%%%%%%%%%%%%%%%
        %%%  Step 1  %%%
        %%%%%%%%%%%%%%%%
        U1 = rand(n,1); %random n-by-1 vector from uniform(0,1)
        R1 = randn(n,1)*sig1; %random n-by-1 vector
        U2 = rand(n,1); %random n-by-1 vector from uniform(0,1)
        R2 = randn(n,1)*sig2; %random
        X1 = [X1;R1(U1 <= kDens(a1,4*bet,R1)./gDens(2*a1+sqrt(8*bet),1/2,R1))];
        X2 = [X2;R2(U2 <= kDens(a2,-4*bet,R2)./gDens(2*a2-4*bet,0,R2))];
        mh = min(length(X1),length(X2));
        X2a = X2(1:mh); X1a=X1(1:mh);
        X2 = X2((X1a.^2+X2a.^2)<1);
        X1 = X1((X1a.^2+X2a.^2)<1);
    end
else % bimodal case
    eta = 1-kappa/2/bet;
    a1 = -2*bet*eta;
    a2 = kappa+2*bet;
    kDens = @(a0,be0,u) (exp(2*(-a0*u.^2-2*be0*u.^4)));
    g1Dens = @(A,B,y0,u) ((exp(-2*(A*(u-y0).^2-B))));
    y0 = sqrt(eta/2);
    B = bet*eta^2/2;
    A = eta*bet/2; 
    sig1 = sqrt(1/4/A);
    sig2 = sqrt(1/4/kappa);
    while(length(X2) < n)
        %%%%%%%%%%%%%%%%
        %%%  Step 1  %%%
        %%%%%%%%%%%%%%%%
        U1 = rand(n,1); %random n-by-1 vector from uniform(0,1)
        R1 = randn(n,1)*sig1+y0; %random n-by-1 vector
        X11 = R1(U1 <= kDens(a1,bet,R1)./g1Dens(A,B,y0,R1));
        U11 = rand(size(X11));
        X11 = X11.*((U11<1/2)-(U11>=1/2));
        X1 = [X1; X11];
        U2 = rand(n,1); %random n-by-1 vector from uniform(0,1)
        if kappa == 0
            R2 =  rand(n,1);
            R21=R2(U2 <= kDens(a2,-bet,R2));
            U21=rand(length(R21),1);
            R21=R21.*((U21<1/2)-(U21>=1/2));
            X2=[X2;R21];  
        else
            R2 = randn(n,1)*sig2; %random
            X2 = [X2;R2(U2 <= kDens(a2,-bet,R2)./g1Dens(a2-2*bet,0,0,R2))]; 
        end
        mh = min(length(X1),length(X2));
        X2a = X2(1:mh); X1a=X1(1:mh);
        X2 = X2((X1a.^2+X2a.^2)<1);
        X1 = X1( (X1a .^2+X2a .^2)<1);
    end
end
X1 = X1(1:n);
X2 = X2(1:n);
CosT = 1-2*(X1.^2+X2.^2);
SinT = sqrt(1-CosT.^2);
U4 = rand(n,1); 
Phi = atan(X2./X1)+(U4>1/2)*pi+Psi;
Y = [cos(Phi).*SinT,sin(Phi).*SinT,CosT];
Y = sign_kappa*Y;

%% rotation (orient the North Pole with mean vector (Mu))
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu);
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end