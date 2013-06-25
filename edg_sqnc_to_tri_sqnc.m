%tri_inds_sqnc = edg_sqnc_to_tri_sqnc(edg_inds_sqnc, edg_vrtx_inds, tri_vrtx_inds)
%
%given a sequence of edges that belong to exactly one triangle, returns an
%array of triangles that own the edges.
%Each triangle in the array appears exactly once.
%It is an error if an edge belongs to more than one triangle.


function tri_inds_sqnc = edg_sqnc_to_tri_sqnc(...
    edg_inds_sqnc, edg_vrtx_inds, tri_vrtx_inds)

tri_inds_sqnc = zeros(size(edg_inds_sqnc));

if numel(tri_inds_sqnc) == 0
    return
end

exstng_tri_ind = edge_ind_to_tri_inds(...
        edg_inds_sqnc(1), edg_vrtx_inds, tri_vrtx_inds);

if numel(exstng_tri_ind) ~= 1
    error(['Each edge in the sequence must '...
           'belong to exactly one triangle']);
end    
    
num_sqnc_tris = 1;
tri_inds_sqnc(num_sqnc_tris) = exstng_tri_ind;


for k=1:numel(edg_inds_sqnc)
    
    exstng_tri_ind = edge_ind_to_tri_inds(...
        edg_inds_sqnc(k), edg_vrtx_inds, tri_vrtx_inds); 
    
    if numel(exstng_tri_ind) ~= 1
        error(['Each edge in the sequence must '...
               'belong to exactly one triangle']);
    end
    
    %add the triangle only if it's not already in the array
    if ~any(exstng_tri_ind == tri_inds_sqnc(1:num_sqnc_tris))
        
        
       num_sqnc_tris = num_sqnc_tris + 1;
       tri_inds_sqnc(num_sqnc_tris) = exstng_tri_ind;
        
        
    end
    
end

tri_inds_sqnc = tri_inds_sqnc(1:num_sqnc_tris); 