function [c ceq grad_c grad_ceq] = sphere_constraint(...
    x, x1x2_mean, cnstrnt_rad)

c      = [];
grad_c = [];

dx = x-x1x2_mean;

dx_norm_sqrd = dx'*dx;

ceq = dx_norm_sqrd - cnstrnt_rad*cnstrnt_rad;

grad_ceq = 2*dx;

end
