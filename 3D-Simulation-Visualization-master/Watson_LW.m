function X = Watson_LW(gamm,n)
% WATSON_LW Generates a random sample of size N from the Watson 
%           distribution (FB4 family)
%
% Y = Watson_LW(gamm,N)
%
% gamm          distribution parameter (-inf < gamma < inf)
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
% "The simulation of spherical distributions in the Fisher-Bingham family"
% by Andrew T.A. Wood (1987)
% Communications in Statistics - Simulation and Computation, 16:3, 885-898,
% DOI: 10.1080/03610918708812624


X=[];

%% Li Wong
if gamm > 0; % bipolar
    ro=4*gamm/(2*gamm+3+sqrt((2*gamm+3)^2-16*gamm));
    r=(3*ro/2/gamm)^3*exp(-3+2*gamm/ro);
    while length(X)< n
        U1= rand(n,1);  %ceil(N/a2)
        U2= rand(n,1);
        S= U1.^2./(1-ro*(1-U1.^2)); % 
        W=gamm*S;
        V=r*U2.^2./((1-ro*S).^3);
        X1 = sqrt(  S(V <= exp(2*W))); %
        V=rand(length(X1),1);
        X1 =X1.*((V<1/2)-(V>=1/2));
        X=[X; X1]; % cos(theta)
    end
else  % gridle
    b=exp(2*gamm)-1;
    
    while length(X)< n
        U1= rand(n,1);  %ceil(N/a2)
        U2= rand(n,1);
        V=log(1+U1*b)/gamm;
        ksi=2*pi*U2;
        c=cos(ksi);
        S1=V.*c.^2;         S2=V-S1;
        Sind=logical((S1<=1).*(S2<=1));
        S1=S1(Sind);          S2=S2(Sind);%[S1 S2]
        d=sqrt(V);
        X1=d.*c; X2=d.*sin(ksi);
        Z=[X1(Sind)';X2(Sind)'];
        
        X=[X; Z(:)];  % cos(theta)
    end
end
X=X(1:n);

