function W2=FB4_Wood(kappa,gamm,n)
% FB4_Wood  Generates a random sample of size n from the Fisher-Bingham_4 distribution
%           as a subfunction that is called by RandFB4gN, RandFB5gN, RandFB6gN
%
% W2 = FB4_Wood(kappa,gamm,N)
%
% kappa         distribution parameter (kappa >= 0)
% gamm          distribution parameter (-inf < gamma < inf)
% N             sample size ( where N is a positive integer)
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik, B.Wainwright
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
% Wood page 888
%%
if kappa==0 && gamm==0;
    U = rand(n,1);
    X = 2*U-1; % cos(theta)
    return
end 
if kappa==0 && gamm~=0;
    X= Watson_LW(gamm,n); % Watson !!!!!!!! which one??
    return
end 
if gamm==0  && kappa ~=0;
    X=vMF_Wood(kappa,n,3);% when gamma==0 FB4 is the von Mises-Fisher distribution
    return
end
%% 
cFB4= @(kappa4,gamm)(integral(@(u)(exp(kappa4*u+gamm*u.^2)),-1,1)); % normalizing constant
W2=[];
%%

if gamm < 0 % if gamma<0, then FB4 is a truncated normal
    
      a0= normcdf((kappa-2*gamm)/sqrt(-gamm*2))-...
        normcdf((kappa+2*gamm)/sqrt(-gamm*2));% acceptance ratio
    a2 = exp(gamm)*cFB4(kappa,gamm)/cFB4(kappa+2*gamm,0); % acceptance ratio
    a1= cFB4(kappa,gamm)/cFB4(kappa,0); % acceptance ratio
    %
    a1p0=kappa*sqrt(-pi/gamm)*exp(-kappa^2/(4*gamm))/(2*sinh(kappa));
    a2p0=a2/a0;
    if max(a1p0,a2p0)<=1  % case 1
        mu1=(-kappa+2*gamm)/sqrt(-2*gamm);
        mu2=(-kappa-2*gamm)/sqrt(-2*gamm);
        q1=1/sqrt(-2*gamm);
        q2=-kappa/2/gamm;
        
        while length(W2)< n
            X2=normrnd(-kappa/2/gamm,sqrt(-1/2/gamm),[ceil(n/a0),1]);
            X2=X2(logical((X2>=-1).*(X2<=1)));
            W2  = [W2;X2 ] ;
        end
       
    elseif   max(1/a1p0,a2p0/a1p0)<=1  % case 2
        r=kappa;
        q1=exp(kappa);
        q2=exp(-kappa);
        while length(W2)< n
            s1= rand(ceil(n/a1),1);%
            u=(log(q1*s1+q2*(1-s1)))/r;
            s2= rand(ceil(n/a1),1);
            u=u(s2<=exp(gamm* u.^2));
            W2  = [W2;u ] ;
        end
    elseif  max(1/a2p0,a1p0/a2p0)<=1 % case 3
        r=kappa+2*gamm;
        q1=exp(r);
        q2=exp(-r);
        while length(W2)< n
            s1= rand(ceil(n/a2),1);
            u=(log(q1*s1+q2*(1-s1)))/r;
            s2= rand(ceil(n/a2),1);
            u= u(s2<=exp(gamm* (1-u).^2));
            W2  = [W2;u ] ;
        end
    end
    
elseif gamm > 0 % corresponds to a positive quadratic term
    W2=[];
    r1=kappa+gamm;
    r2=kappa-gamm;
    m1=exp(r1); m2=exp(-r1); n1=exp(r2); n2=exp(-r2);
    lam=1+exp(-2*gamm);
    p1 = r2*sinh(r1)/(r2*sinh(r1)+r1*sinh(r2));
    
    while length(W2) < n
        s1=rand(ceil(n),1);
        s2=rand(ceil(n),1);
        u=(log(m1*s2+m2*(1-s2)))/r1.*(s1<=p1)+(log(n1*s2+n2*(1-s2)))/r2.*(1-(s1<=p1));
        s3=rand(ceil(n),1);
        x_u=u(s3<=lam*exp(gamm*u.^2)/2./cosh((gamm*u)));
        W2 = [W2;x_u];
    end
end
W2=W2(1:n);