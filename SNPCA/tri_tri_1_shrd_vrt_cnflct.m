function tri_conflict = tri_tri_1_shrd_vrt_cnflct(...
    cand_vert_crds, edge_tri_vert_crds, edge_vert_crds, ...
    tri_vert_crds2, tri_vert_crds3)

%takes two triangles as input that share a single vertex
%the triangles are passed in as coordinates of the vertices in R^2

%The vertex coordinates of the the active edge are stored 
%in edge_tri_vert_crds and edge_vert_crds.
%The coordinates of the proposed triangle's vertices are the coordinates of
%the active edge vertices, and the coordinates of the proposed vertex,
%which are stored in cand_vert_crds.

%The vertex coordinates of the existing triangle are stored in
%edge_tri_vert_crds, tri_vert_crds2, and tri_vert_crds3.

%Loosely speaking, the triagles conflict if their intersection is larger 
%than a single point.

%The algorithm may fail if part of the active edge lies on the interior of
%the existing triangle

%The criteria for determining if the triangles conflict is:

%1)
%the existing triangle lies on the interior of the proposed triangle

%OR

%2)
%the candidate vertex lies in the cone defined by the edges of the existing
%triangle that share a vertex with the active edge, i.e. the cone
%edge: edge_tri_vert_crds - > tri_vert_crds2
%edge: edge_tri_vert_crds - > tri_vert_crds3

%OR

%3)
%The candidate vertex lies in a region defined by the intersection of the
%two sets S1 and S2 (defined below).
%S1 is the cone whose vertex is the nonshared vertex of the active edge,
%and whose boundary contains the two line segments connecting 
%the nonshared vertex of the active edge to two vertices of the existing
%triangle.
%The two vertices of the existing triangle are chosen so that the angle
%between the two line segments is maximized.
%Let v1 and v2 be the vertices that define S1, and let v3 be the remaining
%vertex.
%If v3 is on the same side of the edge connecting v1 and v2 as the
%nonshared vertex, then S2 is the cone with vertex v3, and whose boundary
%contains the line segments connecting v1 and v3, and v2 and v3.
%Otherwise, (if v3 is on the other side), then S2 is the cone with vertex
%at the shared vertex, and whose boundary contains the line segments
%anchored at the cone vertex and oriented in the direction of the line
%segments connecting v2 and v1, and the active edge

%it's the job of the calling function to project the coordinates into a 2D
%plane
assert(numel(cand_vert_crds)==2);
assert(numel(edge_tri_vert_crds)==2);
assert(numel(edge_vert_crds)==2);
assert(numel(tri_vert_crds2)==2);
assert(numel(tri_vert_crds3)==2);


%criterion 2:
trnsltd_cand_vert_crds = [cand_vert_crds(1) - edge_tri_vert_crds(1); ...
                          cand_vert_crds(2) - edge_tri_vert_crds(2)];
    
v1 = [tri_vert_crds2(1) - edge_tri_vert_crds(1); ...
      tri_vert_crds2(2) - edge_tri_vert_crds(2)];

v2 = [tri_vert_crds3(1) - edge_tri_vert_crds(1); ...
      tri_vert_crds3(2) - edge_tri_vert_crds(2)];
  
dot_v1v2     = dot(v1,v2);
norm_v2_sqrd = dot(v2,v2);


v1_perp = v1 - (dot_v1v2/norm_v2_sqrd)*v2;

norm_v1_sqrd = dot(v1,v1);

v2_perp = v2 - (dot_v1v2/norm_v1_sqrd)*v1;


tri_conflict = dot(v1_perp, trnsltd_cand_vert_crds) >= 0;

if tri_conflict
    
    tri_conflict = dot(v2_perp, trnsltd_cand_vert_crds) >= 0;    
end

%just return instead of putting the rest of this function in an if
%statement
if tri_conflict    
    return
end

trnsltd_edge_vert_crds = [edge_vert_crds(1) - edge_tri_vert_crds(1); ...
                          edge_vert_crds(2) - edge_tri_vert_crds(2)];


tri_conflict = dot(v1_perp, trnsltd_edge_vert_crds) >= 0;
                      
if tri_conflict
   
    tri_conflict = dot(v2_perp, trnsltd_edge_vert_crds) >= 0;
    
end

if tri_conflict    
    return
end

                      
%criterion 3
%translate so edge_vert_crds is the origin

cand_cone_bndry1 = [edge_tri_vert_crds(1) - edge_vert_crds(1); ...
                    edge_tri_vert_crds(2) - edge_vert_crds(2)];
cand_cone_bndry1_len = norm(cand_cone_bndry1);                
normlzd_cand_cone_bndry1 = (1/cand_cone_bndry1_len)*cand_cone_bndry1;

cand_cone_bndry2 = [tri_vert_crds2(1) - edge_vert_crds(1); ...
                    tri_vert_crds2(2) - edge_vert_crds(2)];
cand_cone_bndry2_len = norm(cand_cone_bndry2);
normlzd_cand_cone_bndry2 = (1/cand_cone_bndry2_len)*cand_cone_bndry2;


cand_cone_bndry3 = [tri_vert_crds3(1) - edge_vert_crds(1); ...
                    tri_vert_crds3(2) - edge_vert_crds(2)];
cand_cone_bndry3_len = norm(cand_cone_bndry3);
normlzd_cand_cone_bndry3 = (1/cand_cone_bndry3_len)*cand_cone_bndry3;

trnsltd_cand_vert_crds = [cand_vert_crds(1) - edge_vert_crds(1); ...
                          cand_vert_crds(2) - edge_vert_crds(2)];

%set v1, v2, and v3 so that:
%the rays edge_vert_crds, v1 and edge_vert_crds, v2 define the widest cone
%v3 lies inside the cone
dot_prods = [ dot(normlzd_cand_cone_bndry1, normlzd_cand_cone_bndry2) ...
              dot(normlzd_cand_cone_bndry1, normlzd_cand_cone_bndry3) ...
              dot(normlzd_cand_cone_bndry2, normlzd_cand_cone_bndry3) ];
                  
[min_dot_prod min_dot_prod_ind] = min(dot_prods);

switch min_dot_prod_ind
    case 1
        v1 = cand_cone_bndry1;
        v2 = cand_cone_bndry2;
        v3 = cand_cone_bndry3;
    case 2
        v1 = cand_cone_bndry1;
        v2 = cand_cone_bndry3;
        v3 = cand_cone_bndry2;
    case 3
        v1 = cand_cone_bndry2;
        v2 = cand_cone_bndry3;
        v3 = cand_cone_bndry1;
end


%v1_h = plot(...
%    [edge_vert_crds(1) edge_vert_crds(1)+v1(1)], ...
%    [edge_vert_crds(2) edge_vert_crds(2)+v1(2)], 'r-');
%v2_h = plot(...
%    [edge_vert_crds(1) edge_vert_crds(1)+v2(1)], ...
%    [edge_vert_crds(2) edge_vert_crds(2)+v2(2)], 'r-');
%v3_h = plot(...
%    [edge_vert_crds(1) edge_vert_crds(1)+v3(1)], ...
%    [edge_vert_crds(2) edge_vert_crds(2)+v3(2)], 'r-');

%delete([v1_h v1_perp_h]);



%determine if the proposed triangle edge with vertices:
%candidate vertex and the vertex of the active edge that is not shared
%conflicts with the triangle that shares an vertex with the active edge

dot_v1v2     = dot(v1,v2);
norm_v1_sqrd = dot(v1,v1);

v1_perp = v2 - (dot_v1v2/norm_v1_sqrd)*v1;

%v1_perp_h = plot(...
%    [edge_vert_crds(1)+v1(1) edge_vert_crds(1)+v1(1)+v1_perp(1)], ...
%    [edge_vert_crds(2)+v1(2) edge_vert_crds(2)+v1(2)+v1_perp(2)], 'r-')


%is the candidate vertex in the half space defined by the first cone
%boundary?
tri_conflict = dot(v1_perp, trnsltd_cand_vert_crds) >= 0;


if tri_conflict
    %check if the candidate vertex in the half space defined by the second
    %cone boundary?
    norm_v2_sqrd = dot(v2,v2);
    v2_perp      = v1 - (dot_v1v2/norm_v2_sqrd)*v2;

    tri_conflict = dot(v2_perp, trnsltd_cand_vert_crds) >= 0;
    
%    v2_perp_h = plot(...
%    [edge_vert_crds(1)+v2(1) edge_vert_crds(1)+v2(1)+v2_perp(1)], ...
%    [edge_vert_crds(2)+v2(2) edge_vert_crds(2)+v2(2)+v2_perp(2)], 'r-')

end

if tri_conflict

    %determine if triangle vertex on interior of the cone (v3)
    %is inside the triangle with sides v1 and v2
    intr_tri_vert_bary_crds = [v1 v2]\v3;
    
    if sum(intr_tri_vert_bary_crds) >= 1
        num_half_spaces = 3;
    else
        num_half_spaces = 4;
    end
    
    
    %if control reaches here, then the candidate edge under consideration
    %is in the cone defined by v1 and v2
    if num_half_spaces == 3
        trnsltd_cand_vert_crds = ...
            [trnsltd_cand_vert_crds(1) - v1(1); ...
             trnsltd_cand_vert_crds(2) - v1(2)];
                 
         v2v1_diff        = v2 - v1;
         dot_v1_v2v1_diff = dot(v1, v2v1_diff);
         norm_v2v1_diff_sqrd = dot(v2v1_diff, v2v1_diff);
         v2v1_diff_perp   = ...
             v1 - (dot_v1_v2v1_diff/norm_v2v1_diff_sqrd)*v2v1_diff;
         
         tri_conflict = ...
             dot(v2v1_diff_perp, trnsltd_cand_vert_crds) >= 0;
                    
    elseif num_half_spaces == 4
        trnsltd_cand_vert_crds = ...
            [trnsltd_cand_vert_crds(1) - v3(1); ...
             trnsltd_cand_vert_crds(2) - v3(2)];
        
        v1v3_diff               = v1 - v3;
        v2v3_diff               = v2 - v3;

        dot_v1v3_diff_v3    = dot(v1v3_diff, v3);
        norm_v1v3_diff_sqrd = dot(v1v3_diff, v1v3_diff);
        
        v1v3_diff_perp   = v3 - ...
            (dot_v1v3_diff_v3/norm_v1v3_diff_sqrd)*v1v3_diff;
        
        
        tri_conflict = dot(v1v3_diff_perp, trnsltd_cand_vert_crds) >= 0;
        
        if tri_conflict
            
            dot_v2v3_diff_v3    = dot(v2v3_diff, v3);
            norm_v2v3_diff_sqrd = dot(v2v3_diff, v2v3_diff);
            
            v2v3_diff_perp   = v3 - ...
                (dot_v2v3_diff_v3/norm_v2v3_diff_sqrd)*v2v3_diff;
            
            tri_conflict = ...
                dot(v2v3_diff_perp, trnsltd_cand_vert_crds) >= 0;
        end
    end
end


%criterion 1
if ~tri_conflict
%no edge of the existing triangle intersects an edge of the proposed
%triangle id contol reaches here, so there is no conflict, 
%or the existing triangle is a subset of the proposed triangle

%the existing triangle is a subset of the proposed triangle if and only if
%its non shared vertices are in the interior of the proposed triangle
%only need to check one- if control reaches here, no edge of the proposed
%triangle intersected an edge of the existing triangle, so the nonshared
%vertices of the existing triangle are both inside or both outside the
%proposed triangle.
%only test one

trnsltd_vert_crds = [tri_vert_crds2(1) - cand_vert_crds(1); ...
                     tri_vert_crds2(2) - cand_vert_crds(2)];
    
v1 = [edge_tri_vert_crds(1) - cand_vert_crds(1); ...
      edge_tri_vert_crds(2) - cand_vert_crds(2)];

v2 = [edge_vert_crds(1) - cand_vert_crds(1); ...
      edge_vert_crds(2) - cand_vert_crds(2)];

trnsltd_vert_bary_crds = [v1 v2]\trnsltd_vert_crds;

%the test that each barycentric coordinate is >= 0 is unnecessary
tri_conflict = ...
    trnsltd_vert_bary_crds(1)   >= 0 && ...
    trnsltd_vert_bary_crds(2)   >= 0 && ...
    sum(trnsltd_vert_bary_crds) <= 1;
end
