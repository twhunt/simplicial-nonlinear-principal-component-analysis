function [x0_objctv_val, varargout] = x0_objctv(...
    x, dstnc_to_sphr_cntr, sphr_rds, mtrc_dstnc1, mtrc_dstnc2)
%2nd output (if present) is the gradient)

[x_sphr_cntr_dstnc, x_sphr_cntr_dstnc_grad] = dstnc_to_sphr_cntr(x);
[x_mtrc_dstnc1, x_mtrc_dstnc1_grad]         = mtrc_dstnc1(x);
[x_mtrc_dstnc2, x_mtrc_dstnc2_grad]         = mtrc_dstnc2(x);

sphr_cnstrnt        = (x_sphr_cntr_dstnc - sphr_rds)/sphr_rds;
x_mtrc_dstnc_dffrnc = x_mtrc_dstnc1 - x_mtrc_dstnc2;
x_mtrc_dstnc_sum    = x_mtrc_dstnc1 + x_mtrc_dstnc2;
iscls_cnstrnt       = x_mtrc_dstnc_dffrnc/x_mtrc_dstnc_sum;

x0_objctv_val = sphr_cnstrnt^2; 
x0_objctv_val = x0_objctv_val + (iscls_cnstrnt)^2;

if nargout == 2
    
    %contribution from sphere constraint
    varargout{1} = ...
        (sphr_cnstrnt/sphr_rds)...
        *x_sphr_cntr_dstnc_grad;

    %contribution from isosceles constraint
    varargout{1} = varargout{1} ...
        + (2*iscls_cnstrnt/(x_mtrc_dstnc_sum^2))...
        *(x_mtrc_dstnc2*x_mtrc_dstnc1_grad ...
        - x_mtrc_dstnc1*x_mtrc_dstnc2_grad);
    
    varargout{1} = 2*varargout{1};

end