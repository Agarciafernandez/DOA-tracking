function X=vMF_Inv(kappa,n)
% VMFB Generates a random sample of size N from the von Mises-Fisher distribution
%
% W1=vMFb(kappa,N)
%
% kappa         distribution parameter, concentration parameter (kappa >= 0)
% N             sample size ( where N is a positive integer)
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
% "Simulation of the von mises fisher distribution"
% by Andrew T.A. Wood (1994)
% Communications in Statistics - Simulation and Computation, 23:1, 157-164,
% DOI: 10.1080/03610919408813161
% dim of the space / m=3 corresponds to the unit sphere

%% Rubinstein 81, p.39, Fisher 87, p.59

     kappaS=sign(kappa);
     kappa=abs(kappa);
%      kappa1=exp(-2*kappa);
     U = rand(n,1); %random n-by-1 vector from uniform(0,1)
     X=log(2*U*sinh(kappa)+exp(-kappa))/kappa;
%      cos(2*asin(sqrt(-log(U*(1-kappa1)+kappa1)/(2*kappa)))); %

 X=kappaS*X(1:n);
end