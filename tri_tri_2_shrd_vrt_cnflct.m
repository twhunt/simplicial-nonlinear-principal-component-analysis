function overlap = tri_tri_2_shrd_vrt_cnflct(...
    shrd_edge_x, shrd_edge_y, shrd_edge_z, ...
    vert_x, vert_y, vert_z)

%The vertex coordinates of the candidate triangle are 
%(shrd_edge_x(1), shrd_edge_y(1), shrd_edge_z(1)), 
%(shrd_edge_x(2), shrd_edge_y(2), shrd_edge_z(2)), and
%(vert_x, vert_y, vert_z)

%The vertex coordinates of the existing triangle are 
%The vertex coordinates of the candidate triangle are 
%(shrd_edge_x(1), shrd_edge_y(1), shrd_edge_z(1)), 
%(shrd_edge_x(2), shrd_edge_y(2), shrd_edge_z(2)), and the origin


if numel(shrd_edge_x) ~= 2 || numel(shrd_edge_y) ~= 2 ...
        || numel(shrd_edge_z) ~= 2
    error('Input shared edge coordinate vectors should have 2 entries')    
end

if numel(vert_x) ~= 1 || numel(vert_y) ~= 1 || numel(vert_z) ~= 1
    error(['Input candidate vertex coordinates should be passed ' ...
        'in as two scalars'])
end

%project all vertices into the plane containing the candidate triangle
edg_mat      = zeros(3,2);
edg_mat(1,:) = shrd_edge_x(1:end);
edg_mat(2,:) = shrd_edge_y(1:end);
edg_mat(3,:) = shrd_edge_z(1:end);

edg_mat(1,:) = edg_mat(1,:) - vert_x;
edg_mat(2,:) = edg_mat(2,:) - vert_y;
edg_mat(3,:) = edg_mat(3,:) - vert_z;

%after translation, the non-shared vertex of existing triangle has
%coordinates: (-vert_x, -vert_y, -vert_z)
%and the coordinates of the candidate vertex are: (0,0,0)

[Q plnr_shrd_edge] = qr(edg_mat, 0);

planr_nonshrd_exstng_crds = -(Q.'*[vert_x; vert_y; vert_z]);

%translate so the non-shared vertex of the existing triangle falls on the
%origin (in the plane)
plnr_shrd_edge(1,:) = plnr_shrd_edge(1,:) - planr_nonshrd_exstng_crds(1);
plnr_shrd_edge(2,:) = plnr_shrd_edge(2,:) - planr_nonshrd_exstng_crds(2);

% planr_vert_crds = Q.'*[vert_x; vert_y; vert_z];
% planr_vert_crds = -planr_vert_crds;

overlap = plnr_tri_tri_2_shrd_vrt_cnflct(...
    plnr_shrd_edge(1,:), plnr_shrd_edge(2,:),  ...
    - planr_nonshrd_exstng_crds(1), - planr_nonshrd_exstng_crds(2));
