
function Phi=vMF_Circ_Phi(kappaV)
% VMF_CIRC_PHI  Generates a random sample of size n from the von Mises-Fisher distribution, where n is the dimension 
%               of kappaV
%
% W1=vMF_Circ_Phi(kappaV)
%
% kappaV        vector of size n of concentration parameters (kappa >= 0)
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

m=2; % dim of the space / m=3 corresponds to the unit sphere
b = (-2*kappaV +((4*(kappaV.^2)+(m-1)^2).^(1/2)))/(m-1);
x0 = (1-b)./(1+b);
c = kappaV.*x0 + (m-1)*log(1-(x0).^2);
A = (m-1)/2; %shape parameter for Beta distribution
n=length(kappaV);
W1=NaN(n,1);

%%
while(sum(isnan(W1)) >0 )

    %%%%%%%%%%%%%%%%
    %%%  Step 1  %%%
    %%%%%%%%%%%%%%%%
    n1=sum(isnan(W1));
    U = rand(n1,1); %random n-by-1 vector from uniform(0,1)
    Z = betarnd(A,A,n1,1); %random n1-by-1 vector from Beta(A,A)
    W = (1-(1+b(isnan(W1))).*Z)./(1-(1-b(isnan(W1))).*Z);
            
    %%%%%%%%%%%%%%%%
    %%Steps 2 & 3%%%
    %%%%%%%%%%%%%%%%
    Step2 = kappaV(isnan(W1)).*W+(m-1)*log(1-x0(isnan(W1)).*W)-c(isnan(W1));
    W( Step2 <= log(U) )=NaN;
    W1( isnan(W1) )= W; %Good Indices
end
U = randi(4,n,1);
Phi = acos(W1)/2; %+((U<=1/2)-(U>1/2))*pi
Phi = Phi.*((U==1)+(U==3))+ ((U==2)+(U==3))*pi-Phi.*((U==2)+(U==4))+(U==4)*2*pi; % doing 2 pi

