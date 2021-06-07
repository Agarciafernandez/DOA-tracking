function log_norm_const=c_Kent(kappa,beta)
%This function calculates the logarithm of the normalising constant of the Kent function
%based on the method explained in
%Kasarapu "Modelling of directional data using Kent distributions"

%We need to calculate S_1^0
m=0;
delta1=2*pi*sqrt(2/kappa);
term_e=2*beta/kappa;

log_norm_const=Calculate_S(m,delta1,term_e,kappa);
