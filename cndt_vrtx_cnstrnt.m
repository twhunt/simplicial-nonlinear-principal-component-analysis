function [c, ceq, varargout] = cndt_vrtx_cnstrnt(x, ...
    sphr_cntr_dstnc, sphr_rds, mtrc1, mtrc2)

c = [];

[mtrc_val1, mtrc_grad1] = mtrc1(x);
[mtrc_val2, mtrc_grad2] = mtrc2(x);

[sphr_cntr_dstnc_val, sphr_cntr_dstnc_grad] = sphr_cntr_dstnc(x);

mtrc_sum    = mtrc_val1 + mtrc_val2; 
mtrc_dffrnc = mtrc_val1 - mtrc_val2;

%isosceles equality constraint
ceq(2) = mtrc_dffrnc/mtrc_sum;
%sphere constraint
ceq(1) = (sphr_cntr_dstnc_val - sphr_rds)/sphr_rds;

if nargout == 4
    
    %gradient of nonlinear inequality constraints
    varargout{1} = [];
    
    %isosceles constraint gradient
    varargout{2}(:, 2) = mtrc_val2*mtrc_grad1 - mtrc_val1*mtrc_grad2;
    varargout{2}(:, 2) = (2/mtrc_sum^2)*varargout{2}(:, 2);
    
    %sphere constraint gradient
    varargout{2}(:, 1) = (1/sphr_rds)*sphr_cntr_dstnc_grad;
    
end