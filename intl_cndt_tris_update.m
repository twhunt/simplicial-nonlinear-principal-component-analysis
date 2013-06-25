function [...
    new_tri_vrtx_inds, new_edg_vrtx_inds, new_vrtx_inds, ...
    new_vrtx_crdnts, ...
    srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind] ...
    = ...
    intl_cndt_tris_update(...
    intl_tri_is_vbl, intl_edg_srfc_inds, intl_vrtx_inds, ...
    max_vrtx_indx, srfc_crdnts, ...
    srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind)


num_new_tris  = sum(intl_tri_is_vbl);
assert(num_new_tris == 1 || num_new_tris == 2);
crdnts_dmnsn  = size(srfc_crdnts, 1);

if num_new_tris == 1
    
    num_new_vrtxs = 3;
    new_vrtx_inds = max_vrtx_indx + (1:num_new_vrtxs);
        
    new_tri_vrtx_inds = new_vrtx_inds;
    
    new_edg_vrtx_inds = ...
        [ ...
        new_vrtx_inds; ...
        new_vrtx_inds([1 3 2]); ...
        new_vrtx_inds([2 3 1]) ...
        ];
    
else
    %num_new_tris is 2
    
    num_new_vrtxs = 4;
    new_vrtx_inds = max_vrtx_indx + (1:num_new_vrtxs);
    
    new_tri_vrtx_inds = ...
        [...
        new_vrtx_inds(1:3); ...
        new_vrtx_inds([1 2 4]) ...
        ];
    
    new_edg_vrtx_inds = ...
        [...
        new_vrtx_inds(1:3); ...
        new_vrtx_inds([1 3 2]); ...
        new_vrtx_inds([2 3 1]); ...
        new_vrtx_inds([1 4 2]); ...
        new_vrtx_inds([2 4 1])...
        ];
    
end

vrtx_ind_to_srfc_pt_ind(new_vrtx_inds(end)) = 0;
vrtx_ind_to_srfc_pt_ind(new_vrtx_inds(1:2)) = intl_edg_srfc_inds;

new_vrtx_crdnts = zeros(crdnts_dmnsn, num_new_vrtxs);
new_vrtx_crdnts(:, 1:2) = srfc_crdnts(:, intl_edg_srfc_inds);

tmp_indx = 3;
for k=1:2
        
    srfc_pt_ind_to_vrtx_ind(intl_edg_srfc_inds(k)) = new_vrtx_inds(k);
    
    if intl_tri_is_vbl(k)
        
        new_vrtx_crdnts(:, tmp_indx) = ...
            srfc_crdnts(:, intl_vrtx_inds(k));
        
        vrtx_ind_to_srfc_pt_ind(new_vrtx_inds(tmp_indx)) = ...
            intl_vrtx_inds(k);
        
        srfc_pt_ind_to_vrtx_ind(intl_vrtx_inds(k)) = ...
            new_vrtx_inds(tmp_indx);

        
        tmp_indx = tmp_indx + 1;
        
    end
    
end