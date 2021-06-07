function log_ckk=ckk_Kent(kappa,beta,log_ck)
%This function calculates the normalising constant of the Kent function
%based on the method explained in
%Kasarapu "Modelling of directional data using Kent distributions"

%We need to calculate S_1^2
m=2;
delta1=2*pi*sqrt(2/kappa);
term_e=2*beta/kappa;

norm_const=Calculate_S(m,delta1,term_e,kappa);

log_ckk=log_ck+log(exp(norm_const-log_ck)+1/kappa);
