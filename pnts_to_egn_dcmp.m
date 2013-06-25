function [P_dmnt_egnvctrs, P_dmnt_egnvals] = pnts_to_egn_dcmp(...
    cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals, use_sprs_algrthm, eigs_opts)

if use_sprs_algrthm
    %use sparse algorithm
    [P_dmnt_egnvctrs, P_dmnt_egnvals] = pnts_to_egn_dcmp_sprs(...
        cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals, eigs_opts);
else
    %use dense algorithm
    [P_dmnt_egnvctrs, P_dmnt_egnvals] = pnts_to_egn_dcmp_dns(...
        cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals);   
end
