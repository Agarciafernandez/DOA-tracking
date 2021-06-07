function Y = Random_FB6_Wood(kappa,bet,gamm,Mu,Psi,n)
%% 
if kappa==0 && bet==0 && gamm==0
    Y=Random_Uni_Norm(n);
    return
elseif  bet==0 && gamm==0
    Y = Random_vMF_Wood(kappa,Mu,n);
    return
elseif kappa==0 && gamm==0
    Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n);
    return
elseif  kappa==0 && bet==0
    Y = Random_Watsom_Fish(gamm,Mu,n);
    return
elseif  kappa==0 && bet==0
     Y = Random_Watsom_Fish(gamm,Mu,n);
    return
elseif  gamm==0
    Y = Random_FB5_GyT(kappa,bet,Mu,Psi,n);
    return
else  bet==0
    Y = Random_FB4(kappa,gamm,Mu,n);
    return
end 

%% p2 calculation

sign_kappa=sign(kappa)+(kappa==0);
kappa=abs(kappa);
p2= PsiFugv_Wood(gamm-bet,kappa)/( PsiFugv_Wood(gamm-bet,kappa)+exp(-2*bet)*PsiFugv_Wood(gamm+bet,kappa));
% q2=1-p2
%% a4
cFB= integral(@(u)(besseli(0,bet*(1-u.^2)).*...
    exp(kappa*u+gamm*u.^2)),-1,1);
cFB0= integral(@(u)(exp(kappa*u+(gamm-bet)*u.^2)),-1,1);
cFB1= integral(@(u)(exp(kappa*u+(gamm+bet)*u.^2)),-1,1);
a4=2*cFB/(exp(bet)*cFB0+exp(-bet)*cFB1);


%%
X=[];
while length(X) < n
    
    y3_minus  = FB4_Wood(kappa,gamm-bet,ceil(n/a4));
%     FB4RandNb(kappa,gamm-bet,ceil(N/a4));
    
    y3_plus  = FB4_Wood(kappa,gamm+bet,ceil(n/a4));
%     FB4RandNb(kappa,gamm+bet,ceil(N/a4));
    y3_Unif = rand(ceil(n/a4),1);
    y3_mix = y3_minus .* (y3_Unif <= p2)+y3_plus .* (1-(y3_Unif <= p2)) ;
    %target dist. eval at u, g(u;kappa,0,gamma) [x3_mix,x3_plus]
    U=rand(ceil(n/a4),1);
    gFB =besseli(0,bet*(1-U.^2))./cosh(bet*(1-U.^2));
    y3_mix = y3_mix(U<=gFB);
    X = [X;y3_mix];
    
end
X=X(1:n);
%% phi Uniform:calculating the x1 and x2 components;
% U1=rand(N,1);
% U2=rand(N,1);

Phi=vMF_Circ_Phi(bet*(1-X.^2)); %vMF
% if Psi~=0
%     PhivMF=cos(acos(PhivMF)-2*Psi);
% end
cPhi = cos(Phi+Psi);
sPhi = sin(Phi+Psi);
% cPhi=sqrt((1+PhivMF)/2) ;
% cPhi=cPhi.*((U2<1/2)-(U2>=1/2));
% sPhi=sqrt((1-PhivMF)/2) ; 

% sPhi=sPhi.*((U1<1/2)-(U1>=1/2));
%     theta=acos(2*rand(N,1)-1);
sX=sqrt(1-X.^2); %
  Y=sign_kappa*[cPhi.*sX,sPhi.*sX,X];

%% rotation (orient the North Pole with mean vector (mu))
Np=[0,0,1]; % z-axis (North Pole)
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end

