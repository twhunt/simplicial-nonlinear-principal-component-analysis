function tri_conflict = tri_tri_0_shrd_vrt_cnflct(...
    tri1_vert_crds1, tri1_vert_crds2, tri1_vert_crds3, ...
    tri2_vert_crds1, tri2_vert_crds2, tri2_vert_crds3)

%indicates if the two triangles overlap
tri_conflict = false;

%input coordinate vectors may be row or column vectors
%copy them so each column of tri_vert_crds{i} is a coordinate vector
[tri_vert_crds{1:2}] = deal(zeros(2,3));
tri_vert_crds{1}(:,1) = tri1_vert_crds1(1:end);
tri_vert_crds{1}(:,2) = tri1_vert_crds2(1:end);
tri_vert_crds{1}(:,3) = tri1_vert_crds3(1:end);

tri_vert_crds{2}(:,1) = tri2_vert_crds1(1:end);
tri_vert_crds{2}(:,2) = tri2_vert_crds2(1:end);
tri_vert_crds{2}(:,3) = tri2_vert_crds3(1:end);


[tri_edges{1:2}]         = deal(zeros(2,2));
[vert_is_intr{1:2}]      = deal(false(3,1));
[vert_bary_crds{1:2}]    = deal(zeros(3,1));
[trnsltd_vert_crds{1:2}] = deal(zeros(2,1));


for i=1:2
    if i == 1
        T1_ind = 1;
        T2_ind = 2;
    else %i == 2
        T1_ind = 2;
        T2_ind = 1;
    end
    
    %this choice of edges is arbitrary. it may be a good idea to choose
    %the two edges that have the larges angle between them so the
    %following system solving is as stable as possible
    tri_edges{T1_ind} = ...
        [tri_vert_crds{T1_ind}(:,2) - tri_vert_crds{T1_ind}(:,1) ...
         tri_vert_crds{T1_ind}(:,3) - tri_vert_crds{T1_ind}(:,1)];
    
    trnsltd_vert_crds{T2_ind} = ...
        [tri_vert_crds{T2_ind}(:,1) - tri_vert_crds{T1_ind}(:,1) ...
         tri_vert_crds{T2_ind}(:,2) - tri_vert_crds{T1_ind}(:,1) ...
         tri_vert_crds{T2_ind}(:,3) - tri_vert_crds{T1_ind}(:,1)];
    
    %compute barycentric coordinates of each vertex of one triangle with
    %respect to the other
    vert_bary_crds{T2_ind} = tri_edges{T1_ind}\trnsltd_vert_crds{T2_ind};
    
    %determine if each vertex of one triangle is inside the other
    vert_is_intr{T2_ind}(1:end) = ...
        all(vert_bary_crds{T2_ind}  >= 0) & ...
        sum(vert_bary_crds{T2_ind}) <= 1;
    
    if any(vert_is_intr{T2_ind})
        tri_conflict = true;
        return;
    end
end


%if control reaches here, then none of the triangle vertices lie inside the
%other triangle.
%the triangles conflict if and only if an edge of one triangle intersects
%the other
tri_conflict = tri_tri_cnflct_helper(...
    tri_vert_crds{1}(:,1), tri_vert_crds{1}(:,1:2), tri_vert_crds{2});
if ~tri_conflict
    tri_conflict = tri_tri_cnflct_helper(...
        tri_vert_crds{1}(:,2), tri_vert_crds{1}(:,3), tri_vert_crds{2});
end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tri_conflict = tri_tri_cnflct_helper(...
    cmmn_vert_crds, vert_crds, tri_vert_crds)

%all coordinate vectors must be column vectors

tri_conflict = false;

%cone_nomrlzd_edges = zeros(2,3);
%cone_edge_len_sqrd = zeros(3,1);
%cone_edge_perp     = zeros(2,1);
%edges_cnflct       = false(size(vert_crds,2),1);

%translate all coordinates so cmmn_vert_crds is the origin
trnsltd_tri_vert_crds = ...
    [tri_vert_crds(1,:) - cmmn_vert_crds(1);
     tri_vert_crds(2,:) - cmmn_vert_crds(2)];

trnsltd_vert_crds = ...
    [vert_crds(1,:) - cmmn_vert_crds(1);
     vert_crds(2,:) - cmmn_vert_crds(2)];

% cone_edge_len_sqrd = ...
%     [dot(trnsltd_tri_vert_crds(:,1), trnsltd_tri_vert_crds(:,1)); ...
%      dot(trnsltd_tri_vert_crds(:,2), trnsltd_tri_vert_crds(:,2)); ...
%      dot(trnsltd_tri_vert_crds(:,3), trnsltd_tri_vert_crds(:,3))]

 cone_edge_len_sqrd = ...
    [trnsltd_tri_vert_crds(:,1)'* trnsltd_tri_vert_crds(:,1); ...
     trnsltd_tri_vert_crds(:,2)'*trnsltd_tri_vert_crds(:,2); ...
     trnsltd_tri_vert_crds(:,3)'*trnsltd_tri_vert_crds(:,3)];

 
cone_nomrlzd_edges = ...
    [(1/sqrt(cone_edge_len_sqrd(1)))*trnsltd_tri_vert_crds(:,1) ...
     (1/sqrt(cone_edge_len_sqrd(2)))*trnsltd_tri_vert_crds(:,2) ...
     (1/sqrt(cone_edge_len_sqrd(3)))*trnsltd_tri_vert_crds(:,3)];
    
% dot_prods = ...
%     [dot(cone_nomrlzd_edges(:,1), cone_nomrlzd_edges(:,2)) ...
%      dot(cone_nomrlzd_edges(:,1), cone_nomrlzd_edges(:,3)) ...
%      dot(cone_nomrlzd_edges(:,2), cone_nomrlzd_edges(:,3))]

%builtin dot function is slow
dot_prods = ...
    [cone_nomrlzd_edges(:,1)'*cone_nomrlzd_edges(:,2) ...
     cone_nomrlzd_edges(:,1)'*cone_nomrlzd_edges(:,3) ...
     cone_nomrlzd_edges(:,2)'*cone_nomrlzd_edges(:,3)];
 
 
[min_dot_prod, min_dot_prod_ind] = min(dot_prods);

switch min_dot_prod_ind
    case 1
        cone_edge_inds = [1 2 3];
    case 2
        cone_edge_inds = [1 3 2];
    case 3
        cone_edge_inds = [2 3 1];
end

trnsltd_tri_vert_crds = trnsltd_tri_vert_crds(:, cone_edge_inds);
cone_edge_len_sqrd    = cone_edge_len_sqrd(cone_edge_inds);

% edges_dot_prd = ...
%     dot(trnsltd_tri_vert_crds(:,1), trnsltd_tri_vert_crds(:,2))

edges_dot_prd = ...
    trnsltd_tri_vert_crds(:,1)'*trnsltd_tri_vert_crds(:,2);

cone_edge_perp = ...
    trnsltd_tri_vert_crds(:,2) ...
    - (edges_dot_prd/cone_edge_len_sqrd(1))*trnsltd_tri_vert_crds(:,1);


%check if each vertex in trnsltd_vert_crds is in the super cone
edges_cnflct = cone_edge_perp.'*trnsltd_vert_crds >= 0;

if any(edges_cnflct)
    cone_edge_perp = ...
        trnsltd_tri_vert_crds(:,1) ...
        - (edges_dot_prd/cone_edge_len_sqrd(2))*trnsltd_tri_vert_crds(:,2);

    edges_cnflct(edges_cnflct) = ...
        cone_edge_perp.'*trnsltd_vert_crds(:, edges_cnflct) >= 0;
else
    %no vertex was in the super cone, so no edge intersection
    %(tri_conflict initialized to false)
    return;
end

if any(edges_cnflct)

    %trnsltd_vert_crds(:,k) is in the super conex if and only if
    %edges_cnflct(k) is true
    vert_bary_crds =...
        trnsltd_tri_vert_crds(:,1:2)\trnsltd_tri_vert_crds(:,3);
    num_half_spaces_is_3 = sum(vert_bary_crds) >= 1;
    
    if num_half_spaces_is_3
        
        trnsltd_vert_crds(1, edges_cnflct) = ...
            trnsltd_vert_crds(1,edges_cnflct) ...
            - trnsltd_tri_vert_crds(1, edges_cnflct);

        trnsltd_vert_crds(2, edges_cnflct) = ...
            trnsltd_vert_crds(2, edges_cnflct) ...
            - trnsltd_tri_vert_crds(2, edges_cnflct);

        %trnsltd_tri_vert_crds(:,3) holds the minor cone edge that does not
        %coincide with the major cone
        trnsltd_tri_vert_crds(:,3) = ...
            trnsltd_tri_vert_crds(:,2) - trnsltd_tri_vert_crds(:,1);
        
        %edges_dot_prd = ...
        %    dot(trnsltd_tri_vert_crds(:,1), trnsltd_tri_vert_crds(:,3));

        edges_dot_prd = ...
            trnsltd_tri_vert_crds(:,1)'*trnsltd_tri_vert_crds(:,3);

        
        %cone_edge_len_sqrd(3) = ...
        %    dot(trnsltd_tri_vert_crds(:,3), trnsltd_tri_vert_crds(:,3))
        
        cone_edge_len_sqrd(3) = ...
            trnsltd_tri_vert_crds(:,3)'*trnsltd_tri_vert_crds(:,3);

        cone_edge_perp = ...
            trnsltd_tri_vert_crds(:,1) ...
            - (edges_dot_prd/cone_edge_len_sqrd(3)) ...
            *trnsltd_tri_vert_crds(:,3);
        
        edges_cnflct(edges_cnflct) = ...
            cone_edge_perp.'*trnsltd_vert_crds(:, edges_cnflct) >= 0;
        
        tri_conflict = any(edges_cnflct);
    else
        %4 half spaces
        
        %
        trnsltd_vert_crds(1, edges_cnflct) = ...
            trnsltd_vert_crds(1,edges_cnflct) ...
            - trnsltd_tri_vert_crds(1, 3);

        trnsltd_vert_crds(2, edges_cnflct) = ...
            trnsltd_vert_crds(2, edges_cnflct) ...
            - trnsltd_tri_vert_crds(2, 3);
        
        %trnsltd_tri_vert_crds(:,1:2) hold the minor cone edges
        trnsltd_tri_vert_crds(:,1) = ...
            trnsltd_tri_vert_crds(:,1) - trnsltd_tri_vert_crds(:,3);

        trnsltd_tri_vert_crds(:,2) = ...
            trnsltd_tri_vert_crds(:,2) - trnsltd_tri_vert_crds(:,3);
        
        %edges_dot_prd = ...
        %    dot(trnsltd_tri_vert_crds(:,1), trnsltd_tri_vert_crds(:,2))

        edges_dot_prd = ...
            trnsltd_tri_vert_crds(:,1)'*trnsltd_tri_vert_crds(:,2);

        %cone_edge_len_sqrd(1:2) = ...
        %    [dot(trnsltd_tri_vert_crds(:,1), trnsltd_tri_vert_crds(:,1));...
        %     dot(trnsltd_tri_vert_crds(:,2), trnsltd_tri_vert_crds(:,2))]

         cone_edge_len_sqrd(1:2) = ...
            [trnsltd_tri_vert_crds(:,1)'*trnsltd_tri_vert_crds(:,1);...
             trnsltd_tri_vert_crds(:,2)'*trnsltd_tri_vert_crds(:,2)];

        cone_edge_perp = ... 
            trnsltd_tri_vert_crds(:,2) ...
            - (edges_dot_prd/cone_edge_len_sqrd(1))...
            *trnsltd_tri_vert_crds(:,1);
        
        edges_cnflct(edges_cnflct) = ...
            cone_edge_perp.'*trnsltd_vert_crds(:, edges_cnflct) >= 0;
        if any(edges_cnflct)
                    cone_edge_perp = ... 
                        trnsltd_tri_vert_crds(:,1) ...
                        - (edges_dot_prd/cone_edge_len_sqrd(2))...
                        *trnsltd_tri_vert_crds(:,2);
                    edges_cnflct(edges_cnflct) = ...
                        cone_edge_perp.' ...
                        *trnsltd_vert_crds(:, edges_cnflct) >= 0;

                    tri_conflict = any(edges_cnflct);
        end
    end
else
    %(tri_conflict initialized to false)
    return;    
end

end

            




