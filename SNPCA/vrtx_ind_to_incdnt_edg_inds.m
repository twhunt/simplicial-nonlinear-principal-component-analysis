function incdnt_edg_inds = vrtx_ind_to_incdnt_edg_inds(...
    vrtx_ind, edg_inds, edg_vrtx_inds)

%build array of logicals where edg_is_incdnt(k) is true means that
%the edge indexed by edge_inds(k) has the vertex indexed by vrtx_ind as one
%of its vertices
edg_is_incdnt =  any(vrtx_ind == edg_vrtx_inds(edg_inds, 1:2), 2);

incdnt_edg_inds = edg_inds(edg_is_incdnt);