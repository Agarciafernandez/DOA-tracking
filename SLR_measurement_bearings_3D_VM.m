function [A_l,b_l,Omega_l]=SLR_measurement_bearings_3D_VM(meank,Pk,weights,W0,Nx,Nz,kappa,Ap_kappa,x_s)
%We make SLR of conditional moments for Von Mises distribution in 3D for
%one sensor located at x_s

%First we compute the moments of h(x), which projects to the unit sphere
%using sigma-points

%Sigma-point generation
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;

%Transformation through h
dist_sigma_points=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2+(sigma_points(5,:)-x_s(3)).^2);

sigma_points_h=[sigma_points(1,:)-x_s(1);sigma_points(3,:)-x_s(2);sigma_points(5,:)-x_s(3)]./repmat(dist_sigma_points,3,1);

mean_h=sigma_points_h*weights';

mean_hh=zeros(Nz);%E[h*h']
cov_xh=zeros(Nx,Nz);

for j=1:length(weights)
    mean_hh=mean_hh+weights(j)*(sigma_points_h(:,j)*sigma_points_h(:,j)');
    cov_xh=cov_xh+weights(j)*((sigma_points(:,j)-meank)*(sigma_points_h(:,j)-mean_h)');
end

%Now we can compute the moments of the transformed variable (using the VM
%adjustments)
p=3; %p=2 as we are on 2D
z_pred_ukf=Ap_kappa*mean_h;
var_pred_ukf=(Ap_kappa)^2*(mean_hh-mean_h*mean_h')+...
    +Ap_kappa/kappa*eye(Nz)+(1-(Ap_kappa)^2-p*Ap_kappa/kappa)*mean_hh';
var_xz_ukf=Ap_kappa*cov_xh;
%Statistical linearisaltion
A_l=var_xz_ukf'/Pk;
b_l=z_pred_ukf-A_l*meank;
Omega_l=var_pred_ukf-A_l*Pk*A_l';

if(det(Omega_l)<0)
    display('Error')
    pause
end







