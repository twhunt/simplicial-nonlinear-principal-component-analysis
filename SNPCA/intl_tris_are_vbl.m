function [intl_tri_is_vbl] = intl_tris_are_vbl(...
    tri_vrtx_inds, edg_vrtx_inds, vrtx_crdnts, ...
    intl_edg_srfc_inds, intl_vrtx_inds, srfc_crdnts, ...
    non_adj_tri_dist_tol)

%the vertex coordinates of the intial triangles are indexed by 
%intl_edg_srfc_inds and intl_vrtx_inds
%the number of initial triangles is equal to the number of entries in
%intl_vrtx_inds.
%If there are two triangles, intl_edg_srfc_inds is the common edge (i.e.,
%it holds the indices into srfc_crdnts of the two vertices common to both
%triangles)

%Basic sanity check:
%Vertices of shared edge are distinct
%non-shared vertices when there are 2 initial triangles are distinct



num_intl_tris = numel(intl_vrtx_inds);
assert(num_intl_tris <= 2);

intl_tri_is_vbl = false(1, num_intl_tris);
if isempty(intl_tri_is_vbl)
    return
end


if intl_edg_srfc_inds(1)==intl_edg_srfc_inds(2)
    
    warning('Vertices of initial edge are not distinct')
    intl_tri_is_vbl(1:end) = false;
    return
    
end


if num_intl_tris == 2 && intl_vrtx_inds(1) == intl_vrtx_inds(2)
    
    %the non-shared vertex indices are not distinct, i.e. two initial
    %triangles coincide
    intl_tri_is_vbl(2) = false;
    num_intl_tris = 1;
    
end

if num_intl_tris == 1
   
    warning(...
        [...
        'Rejecting potentially viable initial triangle.' ...
        ' intl_tris_are_vbl can only be called with two initial '...
        'triangles due to a bad earlier architecture assumption.']);
    intl_tri_is_vbl = false;
    return
    
    %the bad assumption is that cand_vert_error_tmp expects to be called
    %on a candidate triangle containing an existing front edge and a
    %candidate vertex, where the interior vertex of the front edge is not 
    %the candidate vertex. This isn't the case with a single new triangle.
    
end


%check if candidate vertex coincides with the edge vertex, if so,
%don't bother checking for conflict
for tri_i=1:num_intl_tris
    
    intl_tri_is_vbl(tri_i) = ...
        ~any(intl_vrtx_inds(tri_i) == intl_edg_srfc_inds);
    
end


if all(~intl_tri_is_vbl)
    %all initial triangles are inviable due to non-uniqueness of vertices
    return;
end

% if any(~intl_tri_is_vbl)
%     
%     warning(...
%         'SNPCA:intl_tris_are_vbl_deficient_coding', ...
%         ['Initial candidiate vertex coincides with initial ' ...
%         'edge vertex. Rejecting both initial triangles.']);
%     
%     %mark both triangles are inviable
%     intl_tri_is_vbl(1:end) = false;
%     return;
%     
% end


max_vrtx_ind  = max(edg_vrtx_inds(:));
new_vrtx_inds = max_vrtx_ind+[1 2 3];
tmp_edg_vrtx_inds = ...
    [edg_vrtx_inds; ...
    new_vrtx_inds; ...
    new_vrtx_inds([1 3 2]); ...
    new_vrtx_inds([2 3 1]) ...
    ];

tmp_tri_vrtx_inds = [tri_vrtx_inds; new_vrtx_inds];

is_frnt_edg = ...
    edge_blngs_to_xctly_one_tri(tmp_edg_vrtx_inds, tmp_tri_vrtx_inds);
%tmp_frnt_edg_inds = find(is_frnt_edg);


tmp_num_edgs = size(tmp_edg_vrtx_inds, 1);
actv_frnt_edg_indx = tmp_num_edgs - 2;
tmp_frnt_edg_inds = find(is_frnt_edg);
tmp_vrtx_crdnts = ...
    [vrtx_crdnts, zeros(size(srfc_crdnts,1), 3)];
tmp_vrtx_crdnts(:, [end-2 end-1]) = srfc_crdnts(:, intl_edg_srfc_inds);

%intl_tri_is_vbl = true(1,2);


for tri_i=1:num_intl_tris

    if intl_tri_is_vbl(tri_i)
        %coordinates of interior vertex of new edge
        %add tri_i to the existing triangulation, and treat 
        %srfc_crdnts(:, intl_vrtx_inds(tri_j)) as the coordinates of a
        %newly generated triangle
        %newly created triangle to be tested
        tmp_vrtx_crdnts(:, end) = ...
            srfc_crdnts(:, intl_vrtx_inds(tri_i));
        
        %tri_i = 1 => tri_j = 2
        %tri_i = 2 => tri_j = 1
        tri_j = 1 + mod(tri_i,2);        
        
        
        [intl_tri_is_vbl(tri_i), error_struct] = cand_vert_error_tmp(...
            false, ...
            NaN, ...
            srfc_crdnts(:, intl_vrtx_inds(tri_j)), ...
            actv_frnt_edg_indx, ...
            tmp_frnt_edg_inds, ...
            tmp_tri_vrtx_inds, ...
            tmp_edg_vrtx_inds, ...
            tmp_vrtx_crdnts, ...
            non_adj_tri_dist_tol, ...
            1);
        intl_tri_is_vbl(tri_i) = ~intl_tri_is_vbl(tri_i);
    end
    
end

end