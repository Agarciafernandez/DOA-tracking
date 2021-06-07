function log_cb=cb_Kent(kappa,beta)
%This function calculates the normalising constant of the Kent function
%based on the method explained in
%Kasarapu "Modelling of directional data using Kent distributions"
 
%We need to calculate S_2^0
m=0;
delta1=4*pi/beta*sqrt(2/kappa);
term_e=2*beta/kappa;

log_cb=Calculate_S2(m,delta1,term_e,kappa);
