function X=vMF_Wood(kappa,n,d)
% VMFB Generates a random sample of size N from the von Mises-Fisher distribution
%
% W1=vMF_Wood(kappa,d,n);
%
% kappa         distribution parameter, concentration parameter 
% N             sample size ( where N is a positive integer)
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
%
%
% The methods implemented here are based on
% "Simulation of the von mises fisher distribution"
% by Andrew T.A. Wood (1994)
% Communications in Statistics - Simulation and Computation, 23:1, 157-164,
% DOI: 10.1080/03610919408813161
%% Wood 

% d=3;  dim of the space / m=3 corresponds to the unit sphere

    kappa1=kappa;
    kappa=abs(kappa);
    %%%%%%%%%%%%%%%%
    %%%  Step 0  %%%
    %%%%%%%%%%%%%%%%
    b  = (-2*kappa +sqrt(4*(kappa^2)+(d-1)^2) )/(d-1);
    x0 = (1-b)/(1+b);
    c  = kappa*x0 + (d-1)*log(1-x0^2);
    A  = (d-1)/2; %shape parameter for Beta distribution
    % n  = min(2^12,2*n);

    X=[];
    while(length(X) < n)
        %%%%%%%%%%%%%%%%
        %%%  Step 1  %%%
        %%%%%%%%%%%%%%%%
        U = rand(n,1); %random n-by-1 vector from uniform(0,1)
        Z = betarnd(A,A,n,1); %random n-by-1 vector from Beta(A,A)
        W = (1-(1+b)*Z)./(1-(1-b)*Z);
       %%%%%%%%%%%%%%%%
        %%Steps 2 & 3%%%
        %%%%%%%%%%%%%%%%
        Step2 = kappa*W+(d-1)*log(1-x0*W)-c;
        X=[X;W( Step2 > log(U) )]; %Good Indices
    end
    X=sign(kappa1)*X(1:n);
 