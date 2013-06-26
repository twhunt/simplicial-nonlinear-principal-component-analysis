%adjcnt_frnt_edg_cndt_vrtx_inds = cndt_adjcnt_vrtx_indxs(...
%actv_frnt_edge_ind, frnt_edg_inds, edg_vrtx_inds,
%cndt_vrtx_crdnts ,vrtx_crdnts, adj_vert_max_nudge_dist)
%adjcnt_frnt_edg_cndt_vrtx_inds is a list of vertex indices of vertices
%that:
%1) Belong to a front edge adjacent to the active front edge, but are not
%   shared between the two edges
%2) Whose distance to the candidate vertex is within the tolerance
%   specifed by adj_vert_max_nudge_dist
function adjcnt_frnt_edg_cndt_vrtx_inds = cndt_adjcnt_vrtx_indxs(...
    actv_frnt_edge_ind, frnt_edg_inds, edg_vrtx_inds, ...
    cndt_vrtx_crdnts ,vrtx_crdnts, adj_vert_max_nudge_dist)

%build list of vertices belonging to front edge adjacent to
% active edge, but not shared with the active edge
adjcnt_frnt_edg_nonshrd_vtx_inds = ...
    non_shrd_adj_frnt_edg_vrtx_inds(...
    actv_frnt_edge_ind, frnt_edg_inds, edg_vrtx_inds);

%determine if distance from candidate vertex to vertex that belongs
%to an adjacent front edge and not the front triangle is within the
%tolerance
adjcnt_frnt_edg_cndt_vrtx_inds = [];
num_vrtxs_to_chck = numel(adjcnt_frnt_edg_nonshrd_vtx_inds);
if num_vrtxs_to_chck > 0
    
    adjcnt_frnt_edg_nonshrd_vtx_dstncs = ...
        zeros(1, num_vrtxs_to_chck);
    
    for k=1:num_vrtxs_to_chck
        
        tmp_adj_frnt_edg_vrtx_crdnts = ...
            vrtx_crdnts(:, adjcnt_frnt_edg_nonshrd_vtx_inds(k));
        
        adjcnt_frnt_edg_nonshrd_vtx_dstncs(k) = ...
            norm(...
            tmp_adj_frnt_edg_vrtx_crdnts - cndt_vrtx_crdnts);
    end
    
    adjcnt_frnt_edg_nonshrd_vtx_nearby = ...
        adjcnt_frnt_edg_nonshrd_vtx_dstncs ...
        <= adj_vert_max_nudge_dist;
    
    
    adjcnt_frnt_edg_cndt_vrtx_inds = ...
        adjcnt_frnt_edg_nonshrd_vtx_inds(...
        adjcnt_frnt_edg_nonshrd_vtx_nearby);
end

% num_adj_cndt_vrtxs = numel(adjcnt_frnt_edg_cndt_vrtx_inds);
% num_cand_tris = 1 + num_adj_cndt_vrtxs;
