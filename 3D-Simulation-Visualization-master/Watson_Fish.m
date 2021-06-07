function X = Watson_Fish(gamm,n)
% WATSON_FISH   Generates a random sample of size n from the 
%               Watson-Fisher distribution (in the FB4 family)
%
% X = Watson_Fish(gamm,n)
%
% gamm          distribution parameter (-inf < gamma < inf)
% n             sample size ( where n is a positive integer)
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

%% density

X=[];

if gamm > 0; % bipolar
     C=1/(exp(gamm)-1);
    while length(X)< n
        U1 = rand(n,1); 
        U2 = rand(n,1);
        X0 = log(U1/C+1)/gamm;
        X1 = X0(U2 <= exp(gamm*(X0.^2-X0)));
        V = rand(length(X1),1);
        X1 = X1.*((V<1/2)-(V>=1/2));
        X = [X; X1]; %  cos(theta)
    end
else
     C1 = sqrt(-gamm);
     C2 = atan(C1);
    while length(X)< n
        U1 = rand(n,1);
        U2 = rand(n,1);
        X0 = tan(C2*U1)/C1;
        X1 = X0(U2 <= (1-gamm*X0.^2).*exp(gamm*X0.^2));
        V = rand(length(X1),1);
        X1 = X1.*((V<1/2)-(V>=1/2));
        X = [X; X1]; %  cos(theta)
    end
end
X = X(1:n);
