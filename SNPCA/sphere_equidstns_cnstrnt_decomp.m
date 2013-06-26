function [c ceq grad_c grad_ceq] = sphere_equidstns_cnstrnt_decomp(...
    x, P_dmnt_egnvctrs, P_dmnt_egnvals, biases, ...
    x1x2, x1x2_mean, cnstrnt_rad)



objctv_vals = zeros(1,2);
objctv_grads = zeros(size(P_dmnt_egnvctrs{1},1),1);

for k=1:2
    
    [objctv_vals(k) objctv_grads(:,k)] = dstnc_objctv_decomp(...
        x-x1x2{k}, P_dmnt_egnvctrs{k}, P_dmnt_egnvals{k}, biases(k));

end


dx_mean = x-x1x2_mean;
dx_mean_norm_sqrd = dx_mean.'*dx_mean;

c      = [];
grad_c = [];

ceq(2,1) = dx_mean_norm_sqrd - cnstrnt_rad*cnstrnt_rad;
ceq(1,1) = objctv_vals(1) - objctv_vals(2);

grad_ceq = [(objctv_grads(:,1) - objctv_grads(:,2)) 4*dx_mean];

end