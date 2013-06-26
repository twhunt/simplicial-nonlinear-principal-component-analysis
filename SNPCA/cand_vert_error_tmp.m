function [is_error, error_struct] = cand_vert_error_tmp(...
    is_exstng_vrtx, ...
    vrtx_ind, ...
    cand_vert_coords, ...
    actv_frnt_edg_ind, ...%frnt_cycl_edg_inds, ...
    frnt_edg_inds, ...
    tri_vert_inds, ...
    edg_vert_inds, ...
    vrtx_crdnts, ...
    non_adj_tol, ...
    max_num_nonadj_cnflcts)


%cand_vert_coords:
%coordinates of the candidate vertex

%actv_frnt_edg_ind:
%index of the active front edge

%actv_edg_tri_ind:
%index of the triangle that owns the active front edge
    
%adj_frnt_tri_inds:
%indices of triangles that own a front edge that shares a vertex with the
%active edge

%frnt_tri_inds
%indices of triangles that own a front edge that does not share a vertex
%with the active edge (the interior vertex of the active edge may be shared
%with the active edge)

%tri_vert_inds, edg_vert_inds:
%vertex of indices of the triangles and edges that make up the complex

%tri_coords_x, tri_coords_y, tri_coords_z:
%coordinates of vertices

%non_adj_tol:
%The tolerance that determines if two overlapping triangles conflict

%max_num_nonadj_cnflcts:
%the maximum number of non-adacent triangle conflicts with the candidate
%triangle to check for before exiting the function

%stop_if_adj_cnflct:
%if this flag is set, then check for conflict between the candidate 
%triangle and the non-adjacent front triangles if conflict between the
%candidate triangle and an adjacent front triangle is detected



%The candidate triangle is formed by the active front edge, and the
%candidate vertex

% if actv_frnt_edg_ind == 9
%     disp(actv_frnt_edg_ind);
% end

%REWRITE. old version is at end of this file

%is_error indicates if the proposed vertex should not be placed
%the reason can be extracted from error_struct
%reasons:
%1) proposed triange overlaps with an existing triangle. special case is
%when proposed triangle overlaps existing triangle containing the active
%edge, ie self intersection. there is no fixing this situation
%2) proposed triangle overlaps a triangle that is adjacent to the active
%edge. this may be fixed by moving the proposed vertex to coincide with an
%existing vertex belonging to an edjacent triangle
%3) proposed triangle overlaps with a triangle that is not connected to the
%active edge. fix once no more triangles can be placed
%4) stuff

%condition 3) is defined as overlap between the proposed triangle and the
%front edge, and the two points on the active edge and proposed triangle
%are within a user supplied tolerance
%it is the caller's duty to remove the active edge and adjacent edge from
%the list of edges to check (frnt_edg_inds)

%is_error indicates if the proposed triangle conflicts with an existing
%triangle
%the error struct holds more detailed information about the conflict
%is_error     = false;
error_struct = struct(...
    'self_collision', false, ...
    'actv_edg_is_vbl', true, ...
    'err_id', '', ...
    'adj_tri_collision_inds', [], ...
    'non_adj_tri_collision_inds', []);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START: test for self-collision

%check if candidate vertex is a vertex of the active edge
%uses short circuit AND, vrtx_ind may not be a valid vertex indx if
%is_exstng_vrtx is false
actv_edg_vrt_inds = edg_vert_inds(actv_frnt_edg_ind, :);

error_struct.self_collision = ...
    is_exstng_vrtx && any(vrtx_ind == actv_edg_vrt_inds(1:2));

if error_struct.self_collision
    
    is_error = true;
    
    error_struct.err_id          = 'NLPCA:self_collision';
    error_struct.actv_edg_is_vbl = false;
    %detected self intersection-there is no rectifying this conflict, so
    %don't check for other types of conflict
    return
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START: build projection matrix that projects into the plane containing the
%candidate triangle


%build coordinate vectors of candidate and active edge triangle vertices,
%and translate so the candidate vertex is the origin

%DELETE \/
%edge_vert_coords1 = [...
%    tri_coords_x(actv_edg_vrt_inds(1)); ...
%    tri_coords_y(actv_edg_vrt_inds(1)); ...
%    tri_coords_z(actv_edg_vrt_inds(1))] ...
%    - cand_vert_coords(:);
%DELETE /\
edge_vert_coords1 = ...
    vrtx_crdnts(:, actv_edg_vrt_inds(1)) - cand_vert_coords(:);


%DELETE \/
%edge_vert_coords2 = [...
%    tri_coords_x(actv_edg_vrt_inds(2)); ...
%    tri_coords_y(actv_edg_vrt_inds(2)); ...
%    tri_coords_z(actv_edg_vrt_inds(2))] ...
%     - cand_vert_coords(:);
%DELETE /\
edge_vert_coords2 = ...
    vrtx_crdnts(:, actv_edg_vrt_inds(2)) - cand_vert_coords(:);



[Q plnr_edge_vrtx_coords] = qr([edge_vert_coords1 edge_vert_coords2], 0);
%END: build projection matrix that projects into the plane containing the
%candidate triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%translate and project the active edge interior vertex into the plane of
%the candidate triangle

%DELETE \/
% edge_vert_coords_intr = [...
%     tri_coords_x(actv_edg_vrt_inds(3)); ...
%     tri_coords_y(actv_edg_vrt_inds(3)); ...
%     tri_coords_z(actv_edg_vrt_inds(3))] ...
%     - cand_vert_coords(:);
%DELETE /\
edge_vert_coords_intr = ...
    vrtx_crdnts(:, actv_edg_vrt_inds(3)) - cand_vert_coords(:);


plnr_edge_vert_coords_intr = Q.'*edge_vert_coords_intr;

%the planar candidate vertex has coordinates [0;0]
plnr_cand_tri_vrtx_x = [plnr_edge_vrtx_coords(1, :) 0];
plnr_cand_tri_vrtx_y = [plnr_edge_vrtx_coords(2, :) 0];

plnr_exstng_tri_vrtx_x = ...
    [plnr_edge_vrtx_coords(1, :) plnr_edge_vert_coords_intr(1)];

plnr_exstng_tri_vrtx_y = ...
    [plnr_edge_vrtx_coords(2, :) plnr_edge_vert_coords_intr(2)];

if is_exstng_vrtx && vrtx_ind == actv_edg_vrt_inds(3)
    
    %the candidate triangle and the triangle that owns the active edge
    %coincide
    error_struct.self_collision = true;
    
else
    
    tri1_shrd_vertcs = [1 2];
    tri2_shrd_vertcs = [1 2];

    error_struct.self_collision  = plnr_tri_tri_cnflct(...
        plnr_cand_tri_vrtx_x, plnr_cand_tri_vrtx_y, ...
        plnr_exstng_tri_vrtx_x, plnr_exstng_tri_vrtx_y, ...
        tri1_shrd_vertcs, tri2_shrd_vertcs);

end

is_error = error_struct.self_collision;

if is_error
    
    error_struct.err_id          = 'NLPCA:self_collision';
    error_struct.actv_edg_is_vbl = false;
    %detected self intersection-there is no rectifying this conflict, so
    %don't check for other types of conflict
    return
    
end
%END: test for self-collision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %START: find the active edge in the cell front cycles
% 
% num_frnt_cycls = numel(frnt_cycl_edg_inds);
% fnd_frnt_cycl  = false;
% 
% for k=1:num_frnt_cycls
%     
%     fnd_frnt_cycl = any( actv_frnt_edg_ind ==  frnt_cycl_edg_inds{k});
%     
%     if fnd_frnt_cycl
%         actv_frnt_cycl_ind = k;
%         break;
%     end
% end
% 
% if ~fnd_frnt_cycl
%     
%     error(...
%         'NLPCA:cand_vert_error:actv_edg_not_in_frnt_cycle', ...
%         'Did not find the active edge in any of the front cycles');
%     
% end
% 
% %frnt_cycl_edg_inds{actv_frnt_cycl_ind} contains the active edge
% 
% %END: find the active edge in the cell front cycles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


actv_tri_ind = ...
    edge_ind_to_tri_inds(actv_frnt_edg_ind, edg_vert_inds, tri_vert_inds);

assert(numel(actv_tri_ind)==1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START test for conflict with triangles that own a front edge in
% frnt_edg_inds that shares a vertex with the active edge 
%actv_edg_vrtx_is_frnt_edg_vrtx = false(size(frnt_cycl_vrtx_inds));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START: logical test for edge that shares exactly one vertex with the
% active edge

[adj_frnt_edg_inds1 adj_frnt_edg_inds2] = edg_shrs_1_vertex(...
    actv_frnt_edg_ind, frnt_edg_inds, edg_vert_inds);

adj_frnt_edg_inds = union(adj_frnt_edg_inds1, adj_frnt_edg_inds2);

% actv_edg_vrtx_is_frnt_edg_vrtx(:,1) = ...
%     (actv_edg_vrt_inds(1) == frnt_cycl_vrtx_inds(:,1) ...
%     & actv_edg_vrt_inds(2) ~= frnt_cycl_vrtx_inds(:,2)) ...
%     | ...
%     (actv_edg_vrt_inds(1) == frnt_cycl_vrtx_inds(:,2) ...
%     & actv_edg_vrt_inds(2) ~= frnt_cycl_vrtx_inds(:,1));
% 
% actv_edg_vrtx_is_frnt_edg_vrtx(:,2) = ...
%     (actv_edg_vrt_inds(2) == frnt_cycl_vrtx_inds(:,1) ...
%     & actv_edg_vrt_inds(1) ~= frnt_cycl_vrtx_inds(:,2)) ...
%     | ...
%     (actv_edg_vrt_inds(2) == frnt_cycl_vrtx_inds(:,2) ...
%     & actv_edg_vrt_inds(1) ~= frnt_cycl_vrtx_inds(:,1));

% END: logical test for edge that shares exactly one vertex with the
% active edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_frnt_edgs_adj_to_actv_edge = numel(adj_frnt_edg_inds);

adj_frnt_edg_tri_inds  = NaN(1, num_frnt_edgs_adj_to_actv_edge);

%build list of front edges that are adjacent to the active edge
%and list of triangles that own the adjacent front edge

num_adj_tris = 0;
for k=1:num_frnt_edgs_adj_to_actv_edge

    adj_frnt_edg_tri_ind = edge_ind_to_tri_inds(...
        adj_frnt_edg_inds(k), edg_vert_inds, tri_vert_inds); 
    
    assert(numel(adj_frnt_edg_tri_ind) == 1);
    
    
    if adj_frnt_edg_tri_ind ~= actv_tri_ind
        
        num_adj_tris = num_adj_tris + 1;
        adj_frnt_edg_tri_inds(num_adj_tris) = adj_frnt_edg_tri_ind;
        
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUG plot
%dbg_vrtx_inds = ...
%    edg_vert_inds(adj_frnt_edg_inds(1:num_frnt_edgs_adj_to_actv_edge), 1:2);%
%
%dbg_vrtx_inds = dbg_vrtx_inds([1 3 2 4]);
%
%dbg_hndl = plot3(...
%   tri_coords_x(dbg_vrtx_inds(1:2)), ...
%   tri_coords_y(dbg_vrtx_inds(1:2)), ...
%   tri_coords_z(dbg_vrtx_inds(1:2)), '*m', ...
%   tri_coords_x(dbg_vrtx_inds(3:4)), ...
%   tri_coords_y(dbg_vrtx_inds(3:4)), ...
%   tri_coords_z(dbg_vrtx_inds(3:4)),'*c');
%delete(dbg_hndl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


exstng_tri_vrtx_crds = zeros(size(vrtx_crdnts,1),3);

for k=1:num_adj_tris
    
    
    exstng_tri_vrtx_inds = tri_vert_inds(adj_frnt_edg_tri_inds(k), :);
    
    %DELETE \/
    %exstng_tri_vrtx_crds(1,:) = tri_coords_x(exstng_tri_vrtx_inds);
    %exstng_tri_vrtx_crds(2,:) = tri_coords_y(exstng_tri_vrtx_inds);
    %exstng_tri_vrtx_crds(3,:) = tri_coords_z(exstng_tri_vrtx_inds);
    %DELETE /\
    
    exstng_tri_vrtx_crds(:,1) = vrtx_crdnts(:,exstng_tri_vrtx_inds(1));
    exstng_tri_vrtx_crds(:,2) = vrtx_crdnts(:,exstng_tri_vrtx_inds(2));
    exstng_tri_vrtx_crds(:,3) = vrtx_crdnts(:,exstng_tri_vrtx_inds(3));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEBUG PLOT
    %dbg_h = plot3(...
    %    exstng_tri_vrtx_crds(1,:), ...
    %    exstng_tri_vrtx_crds(2,:), ...
    %    exstng_tri_vrtx_crds(3,:), 'm*');
    %delete(dbg_h);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    exstng_tri_vrtx_crds(:,1) = ...
        exstng_tri_vrtx_crds(:,1) - cand_vert_coords(:);
    
    exstng_tri_vrtx_crds(:,2) = ...
        exstng_tri_vrtx_crds(:,2) - cand_vert_coords(:);

    exstng_tri_vrtx_crds(:,3) = ...
        exstng_tri_vrtx_crds(:,3) - cand_vert_coords(:);

    
    plnr_exstng_tri_vrtx_crds = Q.'*exstng_tri_vrtx_crds;
        
    %the candidate triangle and the existing triangle share two
    %vertices in common, one of which is the candidate vertex
    
    %the candidate triangle and the existing triangle share a
    %a single vertex in common, and it is not  the candidate vertex
    
    actv_edg_vrtx_is_exstng_tri_vrtx_ind = ...
        actv_edg_vrt_inds(1) == exstng_tri_vrtx_inds;
    
    if any(actv_edg_vrtx_is_exstng_tri_vrtx_ind)
        
        assert(nnz(actv_edg_vrtx_is_exstng_tri_vrtx_ind) == 1);
        
        tri1_shrd_vertcs(1) = 1;
        
    else
        
        actv_edg_vrtx_is_exstng_tri_vrtx_ind = ...
            actv_edg_vrt_inds(2) == exstng_tri_vrtx_inds;
        
        if any(actv_edg_vrtx_is_exstng_tri_vrtx_ind)
            
            assert(nnz(actv_edg_vrtx_is_exstng_tri_vrtx_ind) == 1);
            
            tri1_shrd_vertcs(1) = 2;
            
        end
    end
        
    tri2_shrd_vertcs(1) = find(actv_edg_vrtx_is_exstng_tri_vrtx_ind);

    %the candidate vertex is not a vertex of the active edge
    if isempty(vrtx_ind)
        cand_vrtx_shrd = false(1,3);
    else
        cand_vrtx_shrd = vrtx_ind == exstng_tri_vrtx_inds;
    end
    
    if any(cand_vrtx_shrd)

        num_shrd_vrtcs = 2;
        
        tri1_shrd_vertcs(2) = 3;
        tri2_shrd_vertcs(2) = find(cand_vrtx_shrd);
        
    else
        num_shrd_vrtcs = 1;
    end
    

    cnflct =  plnr_tri_tri_cnflct(...
        plnr_cand_tri_vrtx_x, plnr_cand_tri_vrtx_y, ...
        plnr_exstng_tri_vrtx_crds(1,:), ...
        plnr_exstng_tri_vrtx_crds(2,:), ...
        tri1_shrd_vertcs(1:num_shrd_vrtcs), ...
        tri2_shrd_vertcs(1:num_shrd_vrtcs));

    if cnflct
        is_error = true;
        error_struct.adj_tri_collision_inds(end+1) = ...
            adj_frnt_edg_tri_inds(k);
    end
    
end
% END test for conflict with triangles that own a front edge in
% frnt_edg_inds that shares a vertex with the active edge 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if max_num_nonadj_cnflcts == 0
    return
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START test for conflict between candidate triangle and triangles that
% that contain a front edge adjacent to the active edge, but are not the
% triangle that owns the active edge

%remove the active edge, and edges that share a vertex with the active edge
%from consideration

%is_non_adj_frnt_edg_ind = ...
%    actv_frnt_edg_ind ~= frnt_cycl_edg_inds{actv_frnt_cycl_ind};

%the active edge is adjacent to itself, so remove it
is_non_adj_frnt_edg_ind = actv_frnt_edg_ind ~= frnt_edg_inds;

for k=1:numel(adj_frnt_edg_inds)
    is_non_adj_frnt_edg_ind = is_non_adj_frnt_edg_ind ...
        & adj_frnt_edg_inds(k) ~= frnt_edg_inds;
end

non_adj_frnt_edg_inds = frnt_edg_inds(is_non_adj_frnt_edg_ind);

num_non_adj_frnt_edg_inds = numel(non_adj_frnt_edg_inds);

if num_non_adj_frnt_edg_inds == 0
    return;
end

tri1_shrd_vertcs = zeros(1,3);
tri2_shrd_vertcs = zeros(1,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START: build array of front triangles that do not own an edge adjacent to
% the active edge

frnt_tri_ind = edge_ind_to_tri_inds(...
    actv_frnt_edg_ind, edg_vert_inds, tri_vert_inds);

non_adj_frnt_tri_inds = [];

%for frnt_cycl_ind=1:num_frnt_cycls
    
tmp_non_adj_frnt_tri_inds = edg_sqnc_to_tri_sqnc(...
    frnt_edg_inds, ...
    edg_vert_inds, tri_vert_inds);

is_non_adj_frnt_tri = frnt_tri_ind ~= tmp_non_adj_frnt_tri_inds;

for k=1:numel(adj_frnt_edg_inds)
    
    exstng_tri_ind = ...
        edge_ind_to_tri_inds(...
        adj_frnt_edg_inds(k), edg_vert_inds, tri_vert_inds);
    
    assert(numel(exstng_tri_ind) == 1);
    
    is_non_adj_frnt_tri = is_non_adj_frnt_tri ...
        & exstng_tri_ind ~= tmp_non_adj_frnt_tri_inds;
    
end

tmp_non_adj_frnt_tri_inds = ...
    tmp_non_adj_frnt_tri_inds(is_non_adj_frnt_tri);

tmp_non_adj_tris = numel(tmp_non_adj_frnt_tri_inds);

if ~isempty(tmp_non_adj_frnt_tri_inds)
    non_adj_frnt_tri_inds((end+1):(end+tmp_non_adj_tris)) = ...
        tmp_non_adj_frnt_tri_inds;
end

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG
% if any(tmp_non_adj_frnt_tri_inds ~= non_adj_frnt_tri_inds(1:num_nonadj_frnt_tris))
%     disp('uh oh')
% end
% DEBUG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% END: build array of front triangles that do not own an edge adjacent to
% the active edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


actv_edg_mdpt_crdnts = .5*(edge_vert_coords1 + edge_vert_coords2);

actv_edg_mdpt_frnt_tri_vrtx_dstnc(3) = 0;

for k=1:numel(non_adj_frnt_tri_inds)
    %exstng_tri_ind = non_adj_frnt_tri_inds(1:arr_ind)
    
    exstng_tri_ind = non_adj_frnt_tri_inds(k);
    
    assert(numel(exstng_tri_ind) == 1)

    exstng_tri_vrtx_inds = tri_vert_inds(exstng_tri_ind, :);
    
    %DELETE \/
    %exstng_tri_vrtx_crds(1,:) = tri_coords_x(exstng_tri_vrtx_inds);
    %exstng_tri_vrtx_crds(2,:) = tri_coords_y(exstng_tri_vrtx_inds);
    %exstng_tri_vrtx_crds(3,:) = tri_coords_z(exstng_tri_vrtx_inds);
    %DELETE /\

    exstng_tri_vrtx_crds(:,1) = vrtx_crdnts(:, exstng_tri_vrtx_inds(1));
    exstng_tri_vrtx_crds(:,2) = vrtx_crdnts(:, exstng_tri_vrtx_inds(2));
    exstng_tri_vrtx_crds(:,3) = vrtx_crdnts(:, exstng_tri_vrtx_inds(3));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEBUG PLOT
    %%subplot(2,1,1);
    %dbg_hndl=plot3(...
    %    exstng_tri_vrtx_crds(1,:), ...
    %    exstng_tri_vrtx_crds(2,:), ...
    %    exstng_tri_vrtx_crds(3,:), 'm*');
    %delete(dbg_hndl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:3
        exstng_tri_vrtx_crds(:, m) = ...
            exstng_tri_vrtx_crds(:, m) - cand_vert_coords(:);
    
    end
    
    
    %     for m=1:3
    %         actv_edg_mdpt_frnt_tri_vrtx_dstnc(m) = ...
    %             norm( exstng_tri_vrtx_crds(:,m) - actv_edg_mdpt_crdnts );
    %     end
    %     if all(actv_edg_mdpt_frnt_tri_vrtx_dstnc > srch_rds)
    %         continue
    %     end
    
    
    
    plnr_exstng_tri_vrtx_crds = Q.'*exstng_tri_vrtx_crds;
        
    num_shrd_vrtcs = 0;

    if is_exstng_vrtx
        cand_vrtx_is_shrd = vrtx_ind == exstng_tri_vrtx_inds;
    
        if any(cand_vrtx_is_shrd)
            
            assert(nnz(cand_vrtx_is_shrd) == 1)
            
            num_shrd_vrtcs = 1;
            tri1_shrd_vertcs(1) = 3;
            tri2_shrd_vertcs(1) = find(cand_vrtx_is_shrd);
            
        end
    end
    
    for kk=1:2
        
        cand_vrtx_is_shrd = actv_edg_vrt_inds(kk) == exstng_tri_vrtx_inds;
        
        if any(cand_vrtx_is_shrd)
            
            assert(nnz(cand_vrtx_is_shrd) == 1)

            num_shrd_vrtcs = num_shrd_vrtcs + 1;
            
            tri1_shrd_vertcs(num_shrd_vrtcs) = kk;
            tri2_shrd_vertcs(num_shrd_vrtcs) = find(cand_vrtx_is_shrd);
            
        end
    end
    
    cnflct = plnr_tri_tri_cnflct(...
        plnr_cand_tri_vrtx_x, plnr_cand_tri_vrtx_y, ...
        plnr_exstng_tri_vrtx_crds(1,:), plnr_exstng_tri_vrtx_crds(2,:), ...
        tri1_shrd_vertcs(1:num_shrd_vrtcs), ...
        tri2_shrd_vertcs(1:num_shrd_vrtcs));
    
    if cnflct && num_shrd_vrtcs == 0 
        
        %tmp_h2 = plot3(...
        %    cand_vert_coords(1) + [edge_vert_coords1(1) edge_vert_coords2(1)], ...
        %    cand_vert_coords(2) + [edge_vert_coords1(2) edge_vert_coords2(2)], ...
        %    cand_vert_coords(3) + [edge_vert_coords1(3) edge_vert_coords2(3)], 'om')

        
        
        %tmp_h = plot3(...
        %    cand_vert_coords(1) + exstng_tri_vrtx_crds(1, :), ...
        %    cand_vert_coords(2) + exstng_tri_vrtx_crds(2, :), ...
        %    cand_vert_coords(3) + exstng_tri_vrtx_crds(3, :), 'om')

        
        
        %[nrst_crds1 nrst_crds2 min_dstnc] = tri_tri_nrst_pts(...
        %    zeros(size(edge_vert_coords1)), ...
        %    edge_vert_coords1, ...
        %    edge_vert_coords2, ...
        %    exstng_tri_vrtx_crds(:, 1), ...
        %    exstng_tri_vrtx_crds(:, 2), ...
        %    exstng_tri_vrtx_crds(:, 3))

        [nrst_crds1, nrst_crds2, min_dstnc, exit_flag] ...
            = tri_tri_nrst_pts2(...
            zeros(size(edge_vert_coords1)), ...
            edge_vert_coords1, ...
            edge_vert_coords2, ...
            exstng_tri_vrtx_crds(:, 1), ...
            exstng_tri_vrtx_crds(:, 2), ...
            exstng_tri_vrtx_crds(:, 3));

        if exit_flag < 0
           
            warning(...
                ['Could not find minimum distance between ' ...
                'candidate and existing triangle. Probable bug'])
            
        end
        
        cnflct = min_dstnc <= non_adj_tol;

        %         if ~cnflct
        %            disp('why?')
        %            dbg1_hndl=plot3(...
        %                cand_vert_coords(1) + nrst_crds1(1), ...
        %                cand_vert_coords(2) + nrst_crds1(2), ...
        %                cand_vert_coords(3) + nrst_crds1(3), 'bo')
        %            dbg2_hndl=plot3(...
        %                cand_vert_coords(1) + nrst_crds2(1), ...
        %                cand_vert_coords(2) + nrst_crds2(2), ...
        %                cand_vert_coords(3) + nrst_crds2(3), 'bo')
        %
        %            delete(dbg1_hndl);
        %            delete(dbg2_hndl);
        %         end
        %
        %         delete(tmp_h);
        %         delete(tmp_h2);
        
    end
    
    if cnflct
        is_error = true;
        error_struct.non_adj_tri_collision_inds(end+1) = exstng_tri_ind;
        
        if numel(error_struct.non_adj_tri_collision_inds) ...
                == max_num_nonadj_cnflcts
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEBUG PLOT
    %delete(dbg_hndl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% END test for conflict between candidate triangle and triangles that
% that contain a front edge adjacent to the active edge, but are not the
% triangle that owns the active edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [adj_edg_inds1 adj_edg_inds2] = edg_shrs_1_vertex(...
    actv_edg_ind, edg_inds, edg_vrtx_inds)

actv_edg_ind_in_edg_vrtx_inds = find(actv_edg_ind == edg_inds);

actv_edg_vrtx_inds = edg_vrtx_inds(actv_edg_ind, 1:2);

%list of vertices belonging to edges in edg_inds
edg_vrtx_inds_to_check = edg_vrtx_inds(edg_inds, 1:2);

shares_vrtx = ...
    actv_edg_vrtx_inds(1) == edg_vrtx_inds_to_check(:,1) ...
    | actv_edg_vrtx_inds(1) == edg_vrtx_inds_to_check(:,2);

shares_vrtx(actv_edg_ind_in_edg_vrtx_inds) = false;

adj_edg_inds1 = edg_inds(shares_vrtx);


shares_vrtx = ...
    actv_edg_vrtx_inds(2) == edg_vrtx_inds_to_check(:,1) ...
    | actv_edg_vrtx_inds(2) == edg_vrtx_inds_to_check(:,2);

shares_vrtx(actv_edg_ind_in_edg_vrtx_inds) = false;

adj_edg_inds2 = edg_inds(shares_vrtx);

end