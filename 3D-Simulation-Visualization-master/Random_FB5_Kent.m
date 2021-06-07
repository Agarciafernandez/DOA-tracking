function Y = Random_FB5_Kent(kappa,bet,varargin)
% RANDOM_FB5_KENT Generates a random sample of size n from the Kent 
%                 distribution, which is in the Fisher-Bingham_5 family of
%                 distributions.
%
% Y = Random_FB5_Kent(kappa,bet,Mu,Psi,n)
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
% Y = Random_FB5_Kent(kappa,bet,Mu,Psi,n)
% Y = Random_FB5_Kent(kappa,bet,Mu,Psi)
% Y = Random_FB5_Kent(kappa,bet,Mu)
% Y = Random_FB5_Kent(kappa,bet)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
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


%% Density 

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

sign_kappa = sign(kappa)+(kappa==0);
kappa = abs(kappa); % kappa should be non-negative

%%  ERROR Mesage if kappa<2*bet! it does not work but mine
if kappa < 2*bet
     errormessage = 'Error: kappa less than 2*beta';
     error(errormessage)
%     errormessage = 'Error: \kappa less than 2*\beta';
%     error('something:anything',errormessage)  % 
    return
end

% calculating the x3 component (of (x1,x2,x3)) 
% and mixing proportion (p2)
a = 4*(kappa-2*bet);
b = 4*(kappa+2*bet);
gam0 = 8*bet;
lam1 = sqrt(a+2*sqrt(gam0));
lam2 = sqrt(b);
eDens= @(lam0,c0,u) (exp(c0-lam0*u));
c2=b/8/kappa;
kDens= @(a0,lam0,u) (exp(-(a0*u.^2+lam0*u.^4)/2));

%% 
X1 = [];
X2 = [];
while(length(X2) < n)
    %%%%%%%%%%%%%%%%
    %%%  Step 1  %%%
    %%%%%%%%%%%%%%%%
    U1 = rand(n,1); %random n-by-1 vector from uniform(0,1)
    R1 = exprnd(1/lam1,n,1); %random n-by-1 vector
    U2 = rand(n,1); %random n-by-1 vector from uniform(0,1)
    R2 = exprnd(1/lam2,n,1); %random 
    
    %%%%%%%%%%%%%%%%
    %%Steps 2 & 3%%%
    %%%%%%%%%%%%%%%%
    X2 = [X2;R2(U2 <= kDens(b,-gam0,R2)./eDens(lam2,c2,R2))]; 
    X1 = [X1;R1(U1 <= kDens(a,lam1,R1)./eDens(lam1,1,R1) )]; 
    mh = min(length(X1),length(X2));
    X2a = X2(1:mh); X1a=X1(1:mh);
    X2 = X2((X1a.^2+X2a.^2)<1);
    X1 = X1((X1a.^2+X2a.^2)<1);    
end
X1 = X1(1:n);
X2 = X2(1:n);
U3 = rand(n,1); 
X2 = ((U3<1/2)-(U3>=1/2)).*X2;

%%
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