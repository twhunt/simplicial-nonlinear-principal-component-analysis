function y = Q_prdct_dcmp(...
    x, P_dmnt_egnvctrs, P_dmnt_egnvals, bias)

x_dmnt_egnvctrs_crndt_vctr = P_dmnt_egnvctrs'*x;

if isempty(bias)

    P_dmnt_egn_vctr_prdct = ...
        (1./P_dmnt_egnvals).*x_dmnt_egnvctrs_crndt_vctr;    
    
    y = P_dmnt_egnvctrs*P_dmnt_egn_vctr_prdct;    
    
else
    
    P_dmnt_egn_vctr_prdct = ...
        (1./(P_dmnt_egnvals + bias)).*x_dmnt_egnvctrs_crndt_vctr;
    
    
    y = P_dmnt_egnvctrs*P_dmnt_egn_vctr_prdct;

    bias_inv = 1/bias;        

    y = y + bias_inv*(x - P_dmnt_egnvctrs*x_dmnt_egnvctrs_crndt_vctr);
    
end
