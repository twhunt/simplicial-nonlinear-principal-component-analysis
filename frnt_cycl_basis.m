function frnt_cycl_edg_inds = frnt_cycl_basis(edg_vrtx_inds, tri_vrtx_inds)

edge_is_exposed = edge_blngs_to_xctly_one_tri(...
    edg_vrtx_inds, tri_vrtx_inds);

if ~any(edge_is_exposed)
    %there are no front edges
    frnt_cycl_edg_inds = {};
end

adjncy_mtrx = frnt_edgs_to_adjncy_mtrx(...
    edg_vrtx_inds(:,1:2), find(edge_is_exposed));

max_bfs_queue_lngth = nnz(edge_is_exposed) + 1;


spnng_trees = {};
bck_edgs1 = {};
bck_edgs2 = {};
num_spanning_trees = 0;

%while any(adjncy_mtrx(1:end))
%while any(any(adjncy_mtrx))
%use sparsity info to only compare entries that are potentially true
 while any(any(adjncy_mtrx))
    %find a vertex in the front to get the bfs search started
    [strt_vrtx tmp] = ind2sub(size(adjncy_mtrx), find(adjncy_mtrx,1));
    
    num_spanning_trees = num_spanning_trees + 1;
    %generate the spanning tree and set of back edges (breadth first search)
    [spnng_trees{num_spanning_trees}, ...
        bck_edgs1{num_spanning_trees}, ...
        bck_edgs2{num_spanning_trees}] = ...
        brdth_frst_spnng_tree(...
        adjncy_mtrx, strt_vrtx, max_bfs_queue_lngth);
    
    %\/
    %dbg_adjncy_mtrx = adjncy_mtrx;
    %dbg_adjncy_mtrx(spnng_trees{num_spanning_trees}) = false;
    adjncy_mtrx(spnng_trees{num_spanning_trees}) = false;
    %/\
    
    %remove all edges in the spanning tree and back edges
    %adjncy_mtrx = adjncy_mtrx & ~spnng_trees{num_spanning_trees};
    
    %\/
    %all(find(adjncy_mtrx) == find(dbg_adjncy_mtrx))
    %/\
    
    adjncy_mtrx(...
        bck_edgs1{num_spanning_trees}, ...
        bck_edgs2{num_spanning_trees}) = false;
    adjncy_mtrx(...
        bck_edgs2{num_spanning_trees}, ...
        bck_edgs1{num_spanning_trees}) = false;
end

num_basis_cycls = 0;
for tree_ind=1:num_spanning_trees
    
    for back_edge_ind=1:numel(bck_edgs1{tree_ind})
    
        num_basis_cycls = num_basis_cycls + 1;
        
        frnt_cycl_vrt_inds = brdth_frst_tree_srch(...
            spnng_trees{tree_ind}, ...
            bck_edgs1{tree_ind}(back_edge_ind), ...
            bck_edgs2{tree_ind}(back_edge_ind), max_bfs_queue_lngth);
    
        frnt_cycl_edg_inds{num_basis_cycls} = cycl_vrt_inds_to_edg_inds(...
            frnt_cycl_vrt_inds, edg_vrtx_inds(:,1:2));
        
    end
    
end

%one basis cycle for each back edge in the spanning tree
% num_basis_cycls = numel(bck_edgs1);
% 
% frnt_cycl_vrt_inds = cell(1, num_basis_cycls);
% frnt_cycl_edg_inds = cell(1, num_basis_cycls);
% 
% for k =1:num_basis_cycls
%     frnt_cycl_vrt_inds{k} = brdth_frst_tree_srch(...
%         spnng_tree, bck_edgs1(k), bck_edgs2(k), max_bfs_queue_lngth);
%     
%     frnt_cycl_edg_inds{k} = cycl_vrt_inds_to_edg_inds(...
%         frnt_cycl_vrt_inds{k}, edg_vrtx_inds(:,1:2));
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     % DEBUG
%     frnt_cycl_edg_inds{k}
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG
% back_edges_adjacency_mat = adjncy_mtrx & ~spnng_tree;
% [back_edges1 back_edges2] = ...
%     ind2sub(size(back_edges_adjacency_mat), find(triu(back_edges_adjacency_mat)))
% 
% num_basis_cycls = numel(back_edges1);
% 
% frnt_cycl_vrt_inds = cell(1, num_basis_cycls);
% frnt_cycl_edg_inds = cell(1, num_basis_cycls);
% 
% for k =1:num_basis_cycls
%     frnt_cycl_vrt_inds{k} = brdth_frst_tree_srch(...
%         spnng_tree, back_edges1(k), back_edges2(k), max_bfs_queue_lngth);
%     
%     frnt_cycl_edg_inds{k} = cycl_vrt_inds_to_edg_inds(...
%         frnt_cycl_vrt_inds{k}, edg_vrtx_inds(:,1:2));
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     % DEBUG
%     frnt_cycl_edg_inds{k}
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%