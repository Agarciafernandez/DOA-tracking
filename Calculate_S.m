function output=Calculate_S(m,delta1,term_e,kappa)

%See Kasarapu "Modelling of directional data using Kent distributions"

log_bessel=log(besseli(1/2+m,kappa));

if(isinf(log_bessel))
    n=1/2+m; 
    %We apply the expansion for large kappa
    log_bessel=-1/2*log(2*pi*kappa)+kappa+log(1-(4*n^2-1)/(8*kappa));
end


log_f0=log(gamma(1/2))+log_bessel;




sum_tj=1; %First term is f0/f0
sum_tj_1=0; %Sum t_j minus the first term, as indicated in the paper to measure convergence

log_fj_1=log_f0;
j=1;
condition=1;
epsilon=10^(-6); %Epsilon used to stop



while(condition)
    
    ratio_bessel=besseli(2*j+1/2+m,kappa)/besseli(2*j-3/2+m,kappa);
    
    if(isnan(ratio_bessel))
        %We apply the expansion for large kappa
        n=2*j+1/2+m;
        ratio_bessel=1+2*(1-n)/kappa;
        
    end
    
    
    log_fj=log_fj_1+log((j-1/2)/j)+2*log(term_e)+log(ratio_bessel);
    

    
    tj=exp(log_fj-log_f0);
    
    sum_tj=sum_tj+tj;
    
    sum_tj_1=sum_tj_1+tj;
    
    if(tj/sum_tj_1<epsilon)
        break;
    end
    
    log_fj_1=log_fj;
    j=j+1;
end



output=log(delta1)+log_f0+log(sum_tj);

end


