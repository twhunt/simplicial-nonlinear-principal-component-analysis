%tri_inds = edge_ind_to_tri_inds(edge_ind)
%returns list of triangle indices of triangles that contain the edge
%indexed by edge_ind
function tri_inds = edge_ind_to_tri_inds(...
    edge_ind, edge_vert_inds, tri_vert_inds)


vert_inds = edge_vert_inds(edge_ind,1:2);
vert_inds = sort(vert_inds);

tri_vert_inds = sort(tri_vert_inds, 2, 'ascend');

edgei_is_trij     = vert_inds(1) == tri_vert_inds(:,1);
edgei_is_trij_inds = find(edgei_is_trij);

edgem_is_trin      = vert_inds(2) == tri_vert_inds(edgei_is_trij, 2);

tri1_ind           = edgei_is_trij_inds(edgem_is_trin);

edgem_is_trin      = vert_inds(2) == tri_vert_inds(edgei_is_trij, 3);
tri2_ind           = edgei_is_trij_inds(edgem_is_trin);

edgei_is_trij      = vert_inds(1) == tri_vert_inds(:,2);
edgem_is_trin      = vert_inds(2) == tri_vert_inds(edgei_is_trij, 3);
edgei_is_trij_inds = find(edgei_is_trij);
tri3_ind           = edgei_is_trij_inds(edgem_is_trin);

tri_inds = [tri1_ind; tri2_ind; tri3_ind];

%tmp_tri_inds = [tri1_ind; tri2_ind; tri3_ind];

%assert(sort(tmp_tri_inds) == sort(tri_inds));