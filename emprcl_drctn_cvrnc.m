%P = emprcl_drctn_cvrnc(cntr_crds, nghbrhd_pts_crds)
%cntr_crds is a 1D array
%nghbrhd_pts_crds is a matrix of coordinates where each column is a
%coordinate vector
function P = emprcl_drctn_cvrnc(cntr_crds, nghbrhd_pts_crds)

if isempty(cntr_crds)
    error('NLPCA:Empty_empty_crds', 'Empty coordinate vector');
end

num_nghbrhd_pts = size(nghbrhd_pts_crds, 2);

if num_nghbrhd_pts == 0
    error(...
        'NLPCA:no_surf_pts_in_sphr', ...
        'No surface points in search sphere');
end

%pt_dim is the number of entries in the coordinate vector cntr_crds
pt_dim = numel(cntr_crds);
P = zeros(pt_dim);
%v = direction column vector
for k=1:num_nghbrhd_pts

    v           = nghbrhd_pts_crds(:, k) - cntr_crds(:);
    v_norm_sqrd = v.'*v;    
    P           = P + ((1/v_norm_sqrd)*v)*v.';
        
end

P = (1/num_nghbrhd_pts)*P;