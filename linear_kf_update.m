function [mean_ukf_act,var_ukf_act]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R,z_real)

%Linear Kalman filter update with linearisation parameters
%(A_l,b_l,Omega_l)

z_pred_j=A_l*meank+b_l;
Psi_j=Pk*A_l';

S_j=A_l*Pk*A_l'+R+Omega_l;

resta_z_ukf=z_real-z_pred_j;

Gain=Psi_j/S_j;

mean_ukf_act=meank+Gain*resta_z_ukf;
var_ukf_act=Pk-Gain*Psi_j';

var_ukf_act=(var_ukf_act+var_ukf_act')/2;