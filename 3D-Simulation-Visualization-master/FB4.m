function X=FB4(kappa,gamm,n)
% FB4       Generates a random sample of size n from the Fisher-Bingham_4 distribution
%           as a subfunction that is called by Random_FB4,Random_FB5,
%           Random_FB6
%
% W2 = FB4(kappa,gamm,n);
%
% kappa         distribution parameter (-inf < kappa < inf <)
% gamm          distribution parameter (-inf < gamma < inf)
% n             sample size ( where N is a positive integer)
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
%
%
% The methods implemented here are based on
% "The simulation of spherical distributions in the Fisher-Bingham family"
% by Andrew T.A. Wood (1987)
% Communications in Statistics - Simulation and Computation, 16:3, 885-898,
% DOI: 10.1080/03610918708812624

%% 

cFB4= @(kappa4,gamm1)(integral(@(u)(exp(kappa4*u+gamm1*u.^2)),-1,1)); % normalizing constant
% ckg=cFB4(kappa-gamm,0);
%% 
if gamm < 0 % if gamma<0, then FB4 is a truncated normal
    if kappa <= -2*gamm
        X=[];
        while length(X)< n
            a0= normcdf(1,-kappa/gamm/2,sqrt(-1/gamm/2))-...
                normcdf(-1,-kappa/gamm/2,sqrt(-1/gamm/2)); % acceptance ratio
            X=normrnd(-kappa/gamm/2,sqrt(-1/gamm/2),[ceil(n/a0),1]);
            X=X(X>=-1 & X<=1);
        end
    elseif kappa > -2*gamm
        X=[];  
        a2 = exp(gamm)*cFB4(kappa,gamm)/cFB4(kappa+2*gamm,0); % acceptance ratio
        while length(X)< n
            x2 = vMF_Wood(kappa+2*gamm,ceil(n/a2),3);
            x_U= rand(ceil(n/a2),1);
            x2 = x2(x_U<=exp(gamm)*exp(-2*gamm*x2+gamm*x2.^2));
            X  = [X;x2 ] ;
        end
    end
elseif gamm > 0 % corresponds to a positive quadratic term
    p1 = cFB4(kappa+gamm,0)/(cFB4(kappa+gamm,0)+cFB4(kappa-gamm,0)); % mixing proportion
    a3= p1*(1+exp(-2*gamm))*cFB4(kappa,gamm)/cFB4(kappa+gamm,0); % acceptance ratio
    X=[];
    while length(X) < n
        x3_minus  = vMF_Wood(kappa-gamm,ceil(n/a3),3); % marginal vMF distribution
        x3_plus  = vMF_Wood(kappa+gamm,ceil(n/a3),3); % marginal vMF distribution
        x3_U = rand(length(x3_plus),1);
        x3_mix = x3_plus .* (x3_U <= p1)+x3_minus .* (1-(x3_U <= p1));
        %mix corresponds to the mixed density 
        gV=((1+exp(-2*gamm))*exp(gamm*(x3_mix.^2)))./(exp(gamm*x3_mix)+exp(-gamm*x3_mix));
        U=rand(ceil(n/a3),1);
        x3_mix = x3_mix(U<=gV);
        X = [X;x3_mix];
    end
end
X=X(1:n);
end