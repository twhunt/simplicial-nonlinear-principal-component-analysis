function edg_inds = all_frnt_edg_inds(edg_vrtx_inds, tri_vrtx_inds)

edg_blngs_to_1_tri = ...
    edge_blngs_to_xctly_one_tri(edg_vrtx_inds, tri_vrtx_inds) ;

edg_inds = find(edg_blngs_to_1_tri);