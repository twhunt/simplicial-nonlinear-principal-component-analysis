function [intl_srfc_pt_inds, intl_tri_srfc_inds] = gen_init_tris2(...
    srfc_pt_crdnts, intl_srfc_pt_ind1, ...
    edg_lngth, cnstrnt_rds, srch_rds1, srch_rds2, emprcl_drctn_bias, ...
    use_sprs_algrthm, ...
    nnz_egnvals)


data_pt_dim    = size(srfc_pt_crdnts, 1);
intl_pt_crdnts = srfc_pt_crdnts(:, intl_srfc_pt_ind1);

%find surface data point that is edg_lngth away from the initial surface
%data point
%find all points in a box centered at the initial point, then calculate 
%distances of box points to the initial point
%finding points in the box is computationaly cheap
extents = zeros(data_pt_dim, 2);
extents(:, 1) = intl_pt_crdnts - 1*edg_lngth;
extents(:, 2) = intl_pt_crdnts + 1*edg_lngth;

ordnt_in_extnt = false(size(srfc_pt_crdnts));
for k=1:size(ordnt_in_extnt,1)
    
    ordnt_in_extnt(k, :) = ...
        extents(k, 1) <= srfc_pt_crdnts(k, :) ...
        & srfc_pt_crdnts(k, :) <= extents(k, 2);
    
end

srfc_pt_in_box      = all(ordnt_in_extnt, 1);
srfc_pt_in_box_inds = find(srfc_pt_in_box);
num_srfc_pts_in_box = numel(srfc_pt_in_box_inds);

delta_crdnts = zeros(data_pt_dim, num_srfc_pts_in_box);
for k=1:data_pt_dim
    
    delta_crdnts(k, :) = ...
        srfc_pt_crdnts(k, srfc_pt_in_box) - intl_pt_crdnts(k);

end

dstncs = sum(delta_crdnts.^2, 1);
dstncs = abs(dstncs - edg_lngth^2);
[nrst_dstnc nrst_ind] = min(dstncs);
intl_srfc_pt_ind2 = srfc_pt_in_box_inds(nrst_ind);

if intl_srfc_pt_ind1 == intl_srfc_pt_ind2
    
    %only point in neighborhood of surface point indexed by 
    %intl_srfc_pt_ind1 is itself.
    %cannot generate initial edge
    intl_srfc_pt_inds  = [];
    intl_tri_srfc_inds = [];
    return
    
end

%return indices into surface points data of the initial edge
%replaced with intl_srfc_pt_inds
%intl_edg_srfc_inds = [intl_srfc_pt_ind1 intl_srfc_pt_ind2];

edg_vrtx_crdnts1 = srfc_pt_crdnts(:, intl_srfc_pt_ind1);
edg_vrtx_crdnts2 = srfc_pt_crdnts(:, intl_srfc_pt_ind2);

edg_unit_vctr = edg_vrtx_crdnts2 - edg_vrtx_crdnts1;
edg_unit_vctr = (1/norm(edg_unit_vctr))*edg_unit_vctr;

edg_mdpnt_crdnts = .5*(edg_vrtx_crdnts1 + edg_vrtx_crdnts2);


%build empirical direction correlation matrices associated with the edge
%vertices

%DELETE \/
% is_in_ngbrhd = in_nghbrhd(srfc_pt_crdnts, intl_srfc_pt_ind1, srch_rds1);
% is_in_ngbrhd(intl_srfc_pt_ind1) = false;
% num_pts_in_ngbhrd = sum(is_in_ngbrhd);
% P1 = emprcl_drctn_crrltn(...
%     edg_vrtx_crdnts1, srfc_pt_crdnts(:, is_in_ngbrhd));
% P1 = (1/num_pts_in_ngbhrd)*P1;
% P1 = P1 + emprcl_drctn_bias*eye(size(P1));
%DELETE /\

%TEST \/
%compute eigen decomposition of empirical direction covariance matrices
%only compute eigen information associated with dominant eigenvalues
P_factor         = cell(1,2);
P_dmnt_egnvctrs  = cell(1,2);
P_dmnt_egnvals   = cell(1,2);
dbg_dstnc_objctv = cell(1,2);

intl_srfc_pt_inds = [intl_srfc_pt_ind1 intl_srfc_pt_ind2];
srch_rds          = [srch_rds1 srch_rds2];

edg_vrtx_crdnts = cell(1,2);
for k=1:2
    edg_vrtx_crdnts{k} = srfc_pt_crdnts(:, intl_srfc_pt_inds(k));
end

eigs_opts = struct('issym', 1, 'isreal', 1);

num_nghbrhd_pnts = zeros(1,2);

%build functions that compute metrics based on empirical direction
%covariance matrices
dstnc_objctv = cell(1,2);
for k=1:2
    
    is_in_ngbrhd = ...
        in_nghbrhd(srfc_pt_crdnts, intl_srfc_pt_inds(k), srch_rds(k));
    is_in_ngbrhd(intl_srfc_pt_inds(k)) = false;
    
    nghbrhd_is_empty = ~any(is_in_ngbrhd);
    
    if ~nghbrhd_is_empty
        
        [P_dmnt_egnvctrs{k}, P_dmnt_egnvals{k}] = pnts_to_egn_dcmp(...
            edg_vrtx_crdnts{k}(:), ...
            srfc_pt_crdnts(:, is_in_ngbrhd), ...
            nnz_egnvals, ...
            use_sprs_algrthm, ...
            eigs_opts);

        %         % \/ DEFUNCT \/
        %         P_factor{k} = emprcl_drctn_cvrnc_fctr(...
        %             edg_vrtx_crdnts{k}(:), srfc_pt_crdnts(:, is_in_ngbrhd));
        %
        %         num_nghbrhd_pnts(k) = sum(is_in_ngbrhd);
        %
        %         P_mtrx_vctr_prdct = @(x) ...
        %             (1/num_nghbrhd_pnts(k))*(P_factor{k}*(P_factor{k}'*x));
        %
        %
        %         [P_dmnt_egnvctrs{k}, P_dmnt_egnvals{k}] = eigs(...
        %             P_mtrx_vctr_prdct, size(P_factor{k},1), ...
        %             nnz_egnvals, 'la', eigs_opts);
        %         P_dmnt_egnvals{k} = diag(P_dmnt_egnvals{k});
        %         % /\ DEFUNCT /\
        
        dstnc_objctv{k} = @(x) dstnc_objctv_decomp(...
            x-edg_vrtx_crdnts{k}(:), ...
            P_dmnt_egnvctrs{k}, P_dmnt_egnvals{k}, emprcl_drctn_bias);
        
    else
        
        %a search neighborhood around an edge vertex was empty, so cannot
        %generate initial triangles
        %flag error by setting intl_srfc_pt_inds and intl_tri_srfc_inds to
        %empty matrices
        intl_srfc_pt_inds  = [];
        intl_tri_srfc_inds = [];
        return
                
    end

end

x0_objctv_decomp = @(x) objctv1_decomp(...
    x, ...
    edg_vrtx_crdnts, P_dmnt_egnvctrs, P_dmnt_egnvals, ...
    [emprcl_drctn_bias emprcl_drctn_bias]);
%TEST /\

%DELETE \/
% is_in_ngbrhd = in_nghbrhd(srfc_pt_crdnts, intl_srfc_pt_ind2, srch_rds2);
% is_in_ngbrhd(intl_srfc_pt_ind2) = false;
% num_pts_in_ngbhrd = sum(is_in_ngbrhd);
% P2 = emprcl_drctn_crrltn(...
%     edg_vrtx_crdnts2, srfc_pt_crdnts(:, is_in_ngbrhd));
% P2 = (1/num_pts_in_ngbhrd)*P2;
% P2 = P2 + emprcl_drctn_bias*eye(size(P2));
%DELETE /\



%get basis vectors for plane that fits data near the edge's midpoint
%use two significant eigenvectors of empirical direction correlation matrix
%as basis for this plane
%use fmincon to get a point on the sphere that is equidistant from both
%edge vertices when distance is measured by induced norm
%starting point is the point on the constraint sphere that
%intersects the line perpendicular to the intial edge and passing through
%its midpoint



%TEST \/ 
is_in_ngbrhd = ...
    in_nghbrhd_crdnts(srfc_pt_crdnts, edg_mdpnt_crdnts, ...
    .5*norm(edg_vrtx_crdnts1 - edg_vrtx_crdnts2));

num_pts_in_ngbhrd = sum(is_in_ngbrhd);

nghbrhd_is_empty = ~any(is_in_ngbrhd);
if nghbrhd_is_empty

    intl_srfc_pt_inds  = [];
    intl_tri_srfc_inds = [];
    return
    
end

Pmid_factor = emprcl_drctn_cvrnc_fctr(...
    edg_mdpnt_crdnts, srfc_pt_crdnts(:, is_in_ngbrhd));

P_mtrx_vctr_prdct = @(x) ...
    (1/num_pts_in_ngbhrd)*(Pmid_factor*(Pmid_factor'*x));

%copmute two dominant eigendirections 
[Pmid_dmnt_egnvctrs Pmid_dmnt_egnvals] = eigs(...
    P_mtrx_vctr_prdct, size(Pmid_factor,1), ...
    2, 'la', eigs_opts);
%P_dmnt_egnvals{k} = diag(P_dmnt_egnvals{k});


%[egn_vctrs egn_vals] = eig(Pmid);
%[srtd_egn_vals srt_inds] = sort(diag(egn_vals), 'descend');
%orthnrml_bss = egn_vctrs(:, srt_inds(1:2));

%project edge vector into space spanned by Pmid's significant eigenvectors
edg_unit_vctr_prjctn = Pmid_dmnt_egnvctrs.'*edg_unit_vctr;

%build 
edg_unit_vctr_prjctn_orth = edg_unit_vctr_prjctn([2 1]);
edg_unit_vctr_prjctn_orth(1) = - edg_unit_vctr_prjctn_orth(1);

edg_unit_vctr_orth = Pmid_dmnt_egnvctrs*edg_unit_vctr_prjctn_orth;
edg_unit_vctr_orth = (1/norm(edg_unit_vctr_orth))*edg_unit_vctr_orth;
%TEST /\

options = optimset('fmincon');
%options.Algorithm = 'sqp';
options.Algorithm = 'interior-point';
%options.Display   = 'notify';
%options.Display   = 'final';
options.Display    = 'off';
%options.DerivativeCheck = 'on';

%starting points for fmincon
x0      = zeros(data_pt_dim, 2);
x0(:,1) = edg_mdpnt_crdnts + cnstrnt_rds*edg_unit_vctr_orth;
%x0(:,2) = edg_mdpnt_crdnts - cnstrnt_rds*edg_unit_vctr_orth;

%the nonlinear constraints:
%x the distance from x to x1 must equal the distance from x to x2 where
%distances are measured in the inv(P1) and inv(P2) induced metrics
%x must lie on the constraint sphere
%eq_nlin_cnstrnt1 = @(x) abs(dot(x-x1, P1\(x-x1)) - dot(x-x2, P2\(x-x2)));



%x0 must also satisfy the isosceles condition, where
%dot(x-x1, P1\(x-x1)) = dot(x-x2, P2\(x-x2));
%minimize abs(dot(x-x1, P1\(x-x1)) - dot(x-x2, P2\(x-x2)))
%subject to the constraint that x lies on the constraint sphere
%my_nlin_cnstrnt2 = @(x) nlin_cnstrnt2(x, eq_nlin_cnstrnt2);
%my_nlin_cnstrnt2 = @(x) deal(...
%    -1, nlin_cnstrnt2(x, eq_nlin_cnstrnt2))

options.GradObj    = 'on';
options.GradConstr = 'on';

%get initial point on constraint sphere that is equidistant from x1 and x2
%in the induced metrics associated with x1 and x2
%tic

%DELETE \/
%x0_objctv  = @(x) objctv1(x, edg_vrtx_crdnts1, edg_vrtx_crdnts2, P1, P2);
%DELETE /\
x0_cnstrnt = @(x) sphere_constraint(x, edg_mdpnt_crdnts, cnstrnt_rds);

% dbg_x01 = x0(1);
% dbg_x02 = x0(2);

%x01 and x02 lie on the constraint sphere and are equidistant from the edge
%vertices in induced metric
%DELETE \/
% for k=1:2
%     [x0(:,k) min_val exit_flag] = fmincon(...
%         x0_objctv, x0(:,k), [], [], [], [], [], [], x0_cnstrnt, options);
% end
%DELETE /\


%TEST \/
%compute one starting point for fmincon, then reflect across edge midpoint,
%and feed reflection to fmincon to get second starting point

[x0(:,1) min_val exit_flag] = fmincon(...
    x0_objctv_decomp, x0(:,1), ...
    [], [], [], [], [], [], x0_cnstrnt, options);

d = x0(:,1) - edg_mdpnt_crdnts;
x0(:,2) = edg_mdpnt_crdnts - d;

[x0(:,2) min_val exit_flag] = fmincon(...
    x0_objctv_decomp, x0(:,2), ...
    [], [], [], [], [], [], x0_cnstrnt, options);

%for k=1:2
%    [x0(:,k) min_val exit_flag] = fmincon(...
%        x0_objctv_decomp, x0(:,k), ...
%        [], [], [], [], [], [], x0_cnstrnt, options);
%end
%TEST /\

% [x0(2) min_val exit_flag] = fmincon(...
%     x0_objctv, x0(2), [], [], [], [], [], [], x0_cnstrnt, options);


%DELETE \/
% dstnc_objctv = @(x) dstnc_objctv1(x, edg_vrtx_crdnts1, P1);
% lotso_cnstrnts = @(x) sphere_equidstns_cnstrnt(...
%     x, P1, P2, edg_vrtx_crdnts1, edg_vrtx_crdnts2, edg_mdpnt_crdnts, ...
%     cnstrnt_rds);
%DELETE /\

%TEST \/
%dbg_x = zeros(size(P1),1); dbg_x(1) = 1;
% dbg_x = rand(size(P_factor{1},1),1);
% 
% 
% [val1 grad1] = dstnc_objctv(dbg_x);
% [dbg_val1 dbg_grad1] = dbg_dstnc_objctv{1}(dbg_x);
% [dbg_val2 dbg_grad2] = dbg_dstnc_objctv{2}(dbg_x);
% [val1; dbg_val1]
% norm(grad1 - dbg_grad1)/norm(grad1)


lotso_cnstrnts_decomp = @(x) sphere_equidstns_cnstrnt_decomp(...
    x, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    [emprcl_drctn_bias emprcl_drctn_bias], ... 
    edg_vrtx_crdnts, edg_mdpnt_crdnts, ...
    cnstrnt_rds);



% [c ceq grad_c grad_ceq] = lotso_cnstrnts(dbg_x);
% [c_dbg ceq_dbg grad_c_dbg grad_ceq_dbg] = lotso_cnstrnts_decomp(dbg_x);
%TEST /\

%tri_vrtx_crdnts(:,1) and tri_vrtx_crdnts(:,2) are candidate coordinates...
%for new triangle vertices (before being nudged to nearest data point)
%DELETE \/
% tri_vrtx_crdnts = zeros(data_pt_dim, 2);
% for k=1:2
%     [tri_vrtx_crdnts(:,k) min_val exit_flag] = fmincon(...
%         dstnc_objctv, x0(:,k), [], [], [], [], [], [], ...
%         lotso_cnstrnts, options);
% end
%DELETE /\

%TEST \/
tri_vrtx_crdnts = zeros(data_pt_dim, 2);
for k=1:2
    [tri_vrtx_crdnts(:,k) min_val exit_flag] = fmincon(...
        dstnc_objctv{k}, x0(:,k), [], [], [], [], [], [], ...
        lotso_cnstrnts_decomp, options);
end
%TEST /\

%nudge coordinates that solve the minimization problem to nearest data 
%point
intl_tri_srfc_inds = zeros(1,2);
extnts = zeros(data_pt_dim, 2);
for k=1:2
    
    dstnc_to_edg_vrtx = norm(tri_vrtx_crdnts(:,k) - edg_vrtx_crdnts1);
    
    for m=1:data_pt_dim
    
        extnts(m,1) = tri_vrtx_crdnts(m,k) - dstnc_to_edg_vrtx;
        extnts(m,2) = tri_vrtx_crdnts(m,k) + dstnc_to_edg_vrtx;
        
    end
    
    is_in_box   = in_box(srfc_pt_crdnts, extnts);
    num_pts_in_box  = sum(is_in_box);
    
    dlta_crdnts = zeros(data_pt_dim, num_pts_in_box);
    
    for m=1:data_pt_dim
        
        dlta_crdnts(m,:) = ...
            srfc_pt_crdnts(m, is_in_box) - tri_vrtx_crdnts(m,k);
        
    end

    dstnc_to_data_sqrd = sum(dlta_crdnts.^2, 1);
    [min_dstnc_sqrd min_ind] = min(dstnc_to_data_sqrd);

    tmp_inds = find(is_in_box);
    intl_tri_srfc_inds(k) = tmp_inds(min_ind);
    
end

if intl_tri_srfc_inds(1) == intl_tri_srfc_inds(2)
    %two new triangles coincide
    intl_tri_srfc_inds = intl_tri_srfc_inds(1);
end

% hold off
% data_hndl = plot3(...
%     srfc_pt_crdnts(1,:), ...
%     srfc_pt_crdnts(2,:), ...
%     srfc_pt_crdnts(3,:), ...
%     '.g', 'MarkerSize', 4)
% axis equal
% hold on
% edg_hndl = plot3(...
%     [edg_vrtx_crdnts1(1) edg_vrtx_crdnts2(1)], ...
%     [edg_vrtx_crdnts1(2) edg_vrtx_crdnts2(2)], ...
%     [edg_vrtx_crdnts1(3) edg_vrtx_crdnts2(3)], 'k-*');
% 
% mdpt_hndl = plot3(...
%     edg_mdpnt_crdnts(1), ...
%     edg_mdpnt_crdnts(2), ...
%     edg_mdpnt_crdnts(3), 'k*');
% 
% x0_hndl = plot3(x0(1,:), x0(2,:), x0(3,:), 'r*');
% 
% 
% tri_vrtx_crdnts_hndl = plot3(...
%     tri_vrtx_crdnts(1,:), ...
%     tri_vrtx_crdnts(2,:), ...
%     tri_vrtx_crdnts(3,:), 'ob');
% 
% tri_vrtx_crdnts_hndl2 = plot3(...
%     srfc_pt_crdnts(1,intl_tri_srfc_inds), ...
%     srfc_pt_crdnts(2,intl_tri_srfc_inds), ...
%     srfc_pt_crdnts(3,intl_tri_srfc_inds), 'xb');