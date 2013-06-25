function edge_blngs_to_xctly_one_tri_out = ...
    edge_blngs_to_xctly_one_tri(edges, tris)

if isempty(edges) || isempty(tris)
    edge_blngs_to_xctly_one_tri_out = false(0);
    return;
end

num_edges = size(edges,1);
%num_tris  = size(tris,1);

edge_blngs_to_xctly_one_tri_out = false(num_edges,1);

srtd_edgs = sort(edges(:,1:2),2);
srtd_tris = sort(tris,2);

% edgv1_is_triv1 = false(num_tris, 1);
% edgv1_is_triv2 = false(num_tris, 1);
% edgv2_is_triv3 = false(num_tris, 1);
% edgv2_is_triv3 = false(num_tris, 1);

for edge_i = 1:num_edges
    edgv1_is_triv1 = srtd_edgs(edge_i, 1) == srtd_tris(:,1);
    edgv1_is_triv2 = srtd_edgs(edge_i, 1) == srtd_tris(:,2);
    edgv2_is_triv2 = srtd_edgs(edge_i, 2) == srtd_tris(:,2);
    edgv2_is_triv3 = srtd_edgs(edge_i, 2) == srtd_tris(:,3);
    
    belongs_to_tri = ...
        (edgv1_is_triv1 & edgv2_is_triv2) ...
        | (edgv2_is_triv3 & (edgv1_is_triv1 | edgv1_is_triv2));

    edge_blngs_to_xctly_one_tri_out(edge_i) = ...
        numel(find(belongs_to_tri)) == 1;
        
end


% function edge_blngs_to_xctly_one_tri_out = ...
%     edge_blngs_to_xctly_one_tri(edges, tris)
% 
% num_edges = size(edges,1);
% num_tris  = size(tris,1);
% 
% edge_blngs_to_xctly_one_tri_out = false(num_edges,1);
% 
% for edge_i = 1:num_edges
%     num_tri_owners = 0;
%     edge_vinds = sort(edges(edge_i, 1:2));
%     for tri_i = 1:num_tris
%         %sorting vertices of triangle edges simplifies the logical test
%         %whether the current edge belongs to the current triangle
%         tri_edge_vinds1 = sort(tris(tri_i,[1 2]));
%         tri_edge_vinds2 = sort(tris(tri_i,[1 3]));
%         tri_edge_vinds3 = sort(tris(tri_i,[2 3]));
%         
%         %does the current edge belong to the current triangle?
%         belongs_to_tri = ...
%             all(tri_edge_vinds1 == edge_vinds) ...
%             || all(tri_edge_vinds2 == edge_vinds) ...
%             || all(tri_edge_vinds3 == edge_vinds);
%         if belongs_to_tri
%             num_tri_owners = num_tri_owners + 1;
%         end
%         
%         if num_tri_owners == 2
%             %stop looking once ownership by two triangles is established
%             break
%         end
%         
%     end
%     assert(num_tri_owners > 0,'NLPCA:edge_belongs_to_no_tri','Edge does not belong to a triangle');
%     edge_blngs_to_xctly_one_tri_out(edge_i) = num_tri_owners == 1;
% end