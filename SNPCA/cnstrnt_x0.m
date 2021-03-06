function [x0, eqlty_cnstrnt_val, exit_flag] = cnstrnt_x0( ...
    dstnc_objctvs_decomp, ...
    dstnc_to_sphr_cntr, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    actv_edg_vrtx_crdnts, ...
    actv_intr_crdnts, ...
    cnstrnt_sphr_cntr, ...
    cnstrnt_sphr_rds, ...
    emprcl_drctn_crrltn_eval_biases, ...
    x0_optns)


%generate x0 that:
%lies on the constraint sphere,
%is Q-orthogonal to the active edge
%lies on the correct side of the active edge 

actv_edg_vctr = actv_edg_vrtx_crdnts(:,2) - actv_edg_vrtx_crdnts(:,1);

Q_actv_edg_vctr = Q_prdct_dcmp(...
    actv_edg_vctr, ...
    P_dmnt_egnvctrs{1}, ...
    P_dmnt_egnvals{1}, ...
    emprcl_drctn_crrltn_eval_biases(1));

Q_actv_edg_vctr = ...
    (1/norm(Q_actv_edg_vctr))*Q_actv_edg_vctr;

%x0 is Q orhtogonal to the active edge
intr_vrtx_sprhr_cntr = cnstrnt_sphr_cntr(:) - actv_intr_crdnts(:);
x0 = intr_vrtx_sprhr_cntr...
    - (Q_actv_edg_vctr'*intr_vrtx_sprhr_cntr)*Q_actv_edg_vctr;

%scale tmp_x0 so it has length equal to the radius of the
%constraint sphere
x0 = (cnstrnt_sphr_rds/norm(x0))*x0;

if x0'*intr_vrtx_sprhr_cntr < 0
%x0 should point away from interior of the active triangle when
%anchored at the active edge midpoint
    x0 = -x0;
end

x0 = cnstrnt_sphr_cntr + x0;


%\/ x0  by unconstrained minimization \/
crrnt_x0_objctv = @(x) x0_objctv(...
    x, ...
    dstnc_to_sphr_cntr, ...
    cnstrnt_sphr_rds, ...
    dstnc_objctvs_decomp{1}, ...
    dstnc_objctvs_decomp{2});

[x0, eqlty_cnstrnt_val, exit_flag] = ...
    fminunc(crrnt_x0_objctv, x0, x0_optns);