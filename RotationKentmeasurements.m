function z_output=RotationKentmeasurements(z_input,x_k,x_s)


%We obtain the cosines and sines for bearings elevation
dist_2D_sigma_points=sqrt((x_k(1)-x_s(1)).^2+(x_k(3)-x_s(2)).^2);
dist_3D_sigma_points=sqrt(dist_2D_sigma_points.^2+(x_k(5)-x_s(3)).^2);

cos_theta=dist_2D_sigma_points./dist_3D_sigma_points;
sin_theta=(x_k(5)-x_s(3))./dist_3D_sigma_points;

cos_phi=(x_k(1)-x_s(1))./dist_2D_sigma_points;
sin_phi=(x_k(3)-x_s(2))./dist_2D_sigma_points;


Rot=[cos_phi.*cos_theta, -sin_phi,-cos_phi.*sin_theta;...
    sin_phi.*cos_theta,cos_phi,-sin_phi.*sin_theta;...
    sin_theta,zeros(size(cos_phi)),cos_theta];

z_output=Rot*z_input;

