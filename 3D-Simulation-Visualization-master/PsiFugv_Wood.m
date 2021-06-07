    function cKP0= PsiFugv(Psi,kappa)
        if Psi<0
            cKP0=sqrt(-pi/Psi)*exp(-kappa^2/4/Psi)*...
                (normcdf((kappa-2*Psi)/sqrt(-2*Psi))-normcdf((kappa+2*Psi)/sqrt(-2*Psi)));
        elseif Psi>0
            tau1=(kappa+2*Psi)/2/sqrt(Psi);  
            tau2=(kappa-2*Psi)/2/sqrt(Psi);
            int1=integral(@(u)(exp(u.^2)),0,tau1);
            int2=integral(@(u)(exp(u.^2)),0,abs(tau2));
            cKP0=sqrt(1/Psi) *exp(Psi)*(exp(kappa-tau1^2)*int1-sign(tau2)*exp(-kappa-tau2^2)*int2);
        elseif Psi==0
            cKP0=(exp(kappa)-exp(-kappa))/kappa;
        end
