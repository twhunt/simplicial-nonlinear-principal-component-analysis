function [P_dmnt_egnvctrs, P_dmnt_egnvals] = ...
    pnts_to_egn_dcmp_sprs(cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals, eigs_opts)

num_pts_in_ngbrhd = size(lcl_srfc_crdnts, 2);

% try
% 
%     P_factor = emprcl_drctn_cvrnc_fctr(cntr_crdnts, lcl_srfc_crdnts);
% 
% catch excptn
%     
%     if strcmp(excptn.identifier, 'NLPCA:no_surf_pts_in_sphr')
%         rethrow(excptn);
%     end
% end

P_factor = emprcl_drctn_cvrnc_fctr(cntr_crdnts, lcl_srfc_crdnts);

P_mtrx_vctr_prdct = @(x) (1/num_pts_in_ngbrhd)*(P_factor*(P_factor'*x));

[P_dmnt_egnvctrs, tmp_diag_evals] = eigs(...
    P_mtrx_vctr_prdct, size(P_factor,1), ...
    nnz_egnvals, 'la', eigs_opts);

%eigs returns diagonal matrix of eigenvalues, so just keep
%the diagonal
P_dmnt_egnvals = diag(tmp_diag_evals);
