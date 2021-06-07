function [A_l,b_l,Omega_l]=SLR_measurement_Kent(meank,Pk,weights,W0,Nx,mean_Kent_100,cov_Kent_100,x_s)
%We make the SLR of the conditional moments for the Kent distribution
%for one sensor located at x_s

ck_c=mean_Kent_100(1); %ck divided by c

%Sigma-point generation
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;

%We obtain the cosines and sines for azimuth elevation
dist_2D_sigma_points=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2);
dist_3D_sigma_points=sqrt(dist_2D_sigma_points.^2+(sigma_points(5,:)-x_s(3)).^2);

cos_theta=dist_2D_sigma_points./dist_3D_sigma_points;
sin_theta=(sigma_points(5,:)-x_s(3))./dist_3D_sigma_points;

cos_phi=(sigma_points(1,:)-x_s(1))./dist_2D_sigma_points;
sin_phi=(sigma_points(3,:)-x_s(2))./dist_2D_sigma_points;

%We compute the transformed sigma-points, which are 3x3 rotation matrices.
Rot_sp=Calculate_Rot(cos_phi,sin_phi,cos_theta,sin_theta);

E_g=ck_c*(Rot_sp(1:3,:)*weights');
E_R=zeros(3);
C_xg=zeros(Nx,3);
C_g=zeros(3);

for j=1:length(weights)
    Rot_j=reshape(Rot_sp(:,j),3,3);
    
    E_R=E_R+weights(j)*(Rot_j*cov_Kent_100*Rot_j');
    
    gamma1_j_sca=ck_c*Rot_sp(1:3,j); %Scaled gamma1_j
    
    C_xg=C_xg+weights(j)*((sigma_points(:,j)-meank)*(gamma1_j_sca-E_g)');
    C_g=C_g++weights(j)*((gamma1_j_sca-E_g)*(gamma1_j_sca-E_g)');
    
    
end

A_l=C_xg'/Pk;
b_l=E_g-A_l*meank;
Omega_l=C_g+E_R-A_l*Pk*A_l';

if(det(Omega_l)<0)
    display('Error')
    pause
end

end

function Rot=Calculate_Rot(cos_phi,sin_phi,cos_theta,sin_theta)

%We code the rotation matrix into a vector (to speed up computation)
% Rot=[cos_phi.*cos_theta; -sin_phi;-cos_phi.*sin_theta;...
%     sin_phi.*cos_theta;cos_phi;-sin_phi.*sin_theta;...
%     sin_theta;zeros(size(cos_phi));cos_theta];

Rot=[cos_phi.*cos_theta;sin_phi.*cos_theta;sin_theta;...
    -sin_phi;cos_phi;zeros(size(cos_phi));...
    -cos_phi.*sin_theta;-sin_phi.*sin_theta;cos_theta];


end




