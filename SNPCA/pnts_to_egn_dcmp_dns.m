function [P_dmnt_egnvctrs, P_dmnt_egnvals] = pnts_to_egn_dcmp_dns(...
    cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals)

num_pts_in_ngbrhd = size(lcl_srfc_crdnts, 2);

emprcl_drctn_cvrnc_mtrx = ...
    emprcl_drctn_cvrnc(cntr_crdnts, lcl_srfc_crdnts);


% [P_dmnt_egnvctrs, tmp_diag_evals] = ...
%     eig(emprcl_drctn_cvrnc_mtrx, 'nobalance');

%emprcl_drctn_cvrnc_mtrx is symmetric positive semidefinite, but eig may
%return complex eigenvectors when emprcl_drctn_cvrnc_mtrx is singular (it
%happens in practice)
%This can't happen with svd, and right singular vectors are eigenvectors in
%the positive semidefinite case, so use svd.

[P_dmnt_egnvctrs_unused, tmp_diag_evals, P_dmnt_egnvctrs] = ...
    svd(emprcl_drctn_cvrnc_mtrx);

%Put eigenvalues in a 1D vector (instead of the diagonal matrix returned by
%svd)
P_dmnt_egnvals = diag(tmp_diag_evals);
%[P_dmnt_egnvals sort_prmtn] = sort(P_dmnt_egnvals, 'descend');

%only keep the dominant eigenvalues and eigenvectors specified by
%nnz_egnvals
P_dmnt_egnvals  = P_dmnt_egnvals(1:nnz_egnvals);
P_dmnt_egnvctrs = P_dmnt_egnvctrs(:, 1:nnz_egnvals);
% %eig called with nobalance flag produces eigenvectors with non-unit length
% for k=1:size(P_dmnt_egnvctrs)
%     
%     P_dmnt_egnvctrs(:,k) = ...
%         (1/norm(P_dmnt_egnvctrs(:,k)))*P_dmnt_egnvctrs(:,k);
%     
% end

