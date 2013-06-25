%U = emprcl_drctn_crrltn(cntr_crds, nghbrhd_pts_crds)
%cntr_crds is a 1D array
%nghbrhd_pts_crds is a matrix of coordinates where each column is a
%coordinate vector
%P = U*U' where P is the empirical covariance matrix
function U = emprcl_drctn_cvrnc_fctr(cntr_crds, nghbrhd_pts_crds)


if isempty(cntr_crds)
    error('NLPCA:Empty_empty_crds', 'Empty coordinate vector');
end

num_nghbrhd_pts = size(nghbrhd_pts_crds, 2);

if num_nghbrhd_pts == 0
    
    error(...
        'NLPCA:no_surf_pts_in_sphr', ...
        'No surface points in search sphere');
    
end

%U is the matrix of normalized local directions
%subtract off center, then normalize
U = nghbrhd_pts_crds;
for k=1:num_nghbrhd_pts
    
    U(:,k) = U(:,k) - cntr_crds(:);
    nrm_invrs = 1/norm(U(:,k));
    U(:,k) = nrm_invrs*U(:,k);
        
end
