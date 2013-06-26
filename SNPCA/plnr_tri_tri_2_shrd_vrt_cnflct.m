function overlap = plnr_tri_tri_2_shrd_vrt_cnflct(...
    shrd_edge_x, shrd_edge_y, ...
    vert_x, vert_y)

%The vertex coordinates of the candidate triangle are 
%(shrd_edge_x(1), shrd_edge_y(1)), 
%(shrd_edge_x(2), shrd_edge_y(2)), and
%(vert_x, vert_y)

%The vertex coordinates of the existing triangle are 
%The vertex coordinates of the candidate triangle are 
%(shrd_edge_x(1), shrd_edge_y(1)), 
%(shrd_edge_x(2), shrd_edge_y(2)), and the origin

if numel(shrd_edge_x) ~= 2 || numel(shrd_edge_y) ~= 2
    error('Input shared edge coordinate vectors should have 2 entries')    
end

if numel(vert_x) ~= 1 || numel(vert_y) ~= 1
    error(['Input candidate vertex coordinates should be passed ' ...
        'in as two scalars'])
end

%Ensure that coordinate vectors are column vectors
edge_mat      = zeros(2);

edge_mat(1,:) = shrd_edge_x(1:end);
edge_mat(2,:) = shrd_edge_y(1:end);

vrt_crds_cvec = [vert_x; vert_y];

bary_crds = edge_mat\vrt_crds_cvec;

overlap = sum(bary_crds) <= 1;

% return
% 
% [Q R] = qr([tri_edge_x(1) tri_edge_x(2); ...
%             tri_edge_y(1) tri_edge_y(2); ...
%             tri_edge_z(1) tri_edge_z(2)],0);
% 
% plane_vert_crds = Q.'*[vert_x; vert_y; vert_z];
% 
% bary_crds = R\plane_vert_crds;
% 
% overlap = sum(bary_crds) <= 1;