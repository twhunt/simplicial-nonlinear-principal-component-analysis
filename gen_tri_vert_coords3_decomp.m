function [new_vrtx_crdnts, min_val, exit_flag] ...
    = ...
    gen_tri_vert_coords3_decomp(...
    edge_vrtx_crdnts, ...
    intr_vrtx_crdnts, ...
    cnstrnt_sphr_cntr_crdnts, ...
    cnstrnt_sphr_rds, ...
    x0, ...
    sphr_cntr_dstnc, ...
    dstnc_mtrcs, ...
    cndt_vrtx_optns)


%displacement vectors from interior vertex to edge vertices
edge_vrtx_vctrs = [...
    edge_vrtx_crdnts(:,1) - intr_vrtx_crdnts, ...
    edge_vrtx_crdnts(:,2) - intr_vrtx_crdnts];

[edge_vrtx_vctrs_Q, R] = qr(edge_vrtx_vctrs, 0);

edg_vrtx_diff = edge_vrtx_crdnts(:,2) - edge_vrtx_crdnts(:,1);
[edg_vrtx_diff_Q, R] = qr(edg_vrtx_diff);

%the first column of edg_vrtx_diff_Q is aligned with the active front edge,
%the other columns are orthogonal to the edge
num_edg_orthnrml_vctrs = size(edg_vrtx_diff_Q,2) - 1;

%edg_vrtx_diff_Q(:, 2:end) is an orhtonormal basis for the orthogonal
%complement of span(edg_vrtx_diff)
intrsctn_span  = [edg_vrtx_diff_Q(:, 2:end) edge_vrtx_vctrs];

intrsctn_null_bss = null(intrsctn_span);

actv_edg_orthgnl_cnstrnt = ...
    edg_vrtx_diff_Q(:, 2:end)...
    *intrsctn_null_bss(1:num_edg_orthnrml_vctrs, :);

%orthonormal basis for linear constraint matrix, ie matrix that forces
%minimizer to fall on the correct side of the active edge
%make basis for intersection orthonormal (for x0)
[actv_edg_orthgnl_cnstrnt, R] = qr(actv_edg_orthgnl_cnstrnt, 0);

%flip sign on basis vector if it's pointing in the wrong direction (wrong
%being that fmincon interprets it as pointing inside the triangle)
drctn_flip = ...
    actv_edg_orthgnl_cnstrnt.'...
    *(cnstrnt_sphr_cntr_crdnts(:) - intr_vrtx_crdnts(:)) >= 0;

actv_edg_orthgnl_cnstrnt(:, drctn_flip) = ...
    -actv_edg_orthgnl_cnstrnt(:, drctn_flip);

lnr_inqlty_rhs = ...
    actv_edg_orthgnl_cnstrnt.'*cnstrnt_sphr_cntr_crdnts;
lnr_inqlty_rhs = lnr_inqlty_rhs(:);

%options = optimset('fmincon');
%options.MaxIter = 4096;
%options.Algorithm = 'sqp';
cndt_vrtx_optns.Algorithm = 'interior-point';
%options.Display   = 'notify';
%options.Display   = 'final';
cndt_vrtx_optns.Display    = 'off';


cndt_vrtx_optns.GradObj    = 'on';
cndt_vrtx_optns.GradConstr = 'on';

cndt_vrtx_optns.Hessian = 'fin-diff-grads';
cndt_vrtx_optns.SubproblemAlgorithm = 'cg';

cndt_vrtx_optns.MaxIter = 2048;
cndt_vrtx_optns.Diagnostics = 'off';
%options.DerivativeCheck = 'on';

%get initial point on constraint sphere that is equidistant from x1 and x2
%in the induced metrics associated with x1 and x2
%tic

% actv_edg_orth_bss_crdnts = rand(size(actv_edg_orthgnl_cnstrnt, 2), 1);
% actv_edg_orth_bss_crdnts = ...
%     (cnstrnt_sphr_rds/norm(actv_edg_orth_bss_crdnts))...
%     *actv_edg_orth_bss_crdnts;


% dbg_edg_h = plot3(...
%     edge_vrtx_crdnts(1 ,:), ...
%     edge_vrtx_crdnts(2 ,:), ...
%     edge_vrtx_crdnts(3 ,:), 'r*');
% 
% dbg_intr_h = plot3(...
%     intr_vrtx_crdnts(1), ...
%     intr_vrtx_crdnts(2), ...
%     intr_vrtx_crdnts(3), 'b*');
% 
% dbg_cntr_h = plot3(...
%     cnstrnt_sphr_cntr_crdnts(1), ...
%     cnstrnt_sphr_cntr_crdnts(2), ...
%     cnstrnt_sphr_cntr_crdnts(3), 'rx')
% 
% delete(dbg_edg_h);
% delete(dbg_intr_h);
% delete(dbg_cntr_h);
% 
% 
% dbg_cnstrnt_dsplcmnt_vctr = ...
%     [cnstrnt_sphr_cntr_crdnts x0];
% 
% dbg_cnstrnt_h = plot3(...
%     dbg_cnstrnt_dsplcmnt_vctr(1, :), ...
%     dbg_cnstrnt_dsplcmnt_vctr(2, :), ...
%     dbg_cnstrnt_dsplcmnt_vctr(3, :), 'r-')
% 
% dbg_x0 = plot3(x0(1), x0(2), x0(3), 'ro')
% 
% delete(dbg_cnstrnt_h);
% delete(dbg_x0);


%TEST \/

%disp(num2str(x0_objctv_decomp(x0)))


%dbg_h1 = plot3(x0(1), x0(2), x0(3), 'rx')


sphr_isos_cnstrnt = @(x) cndt_vrtx_cnstrnt(x, ...
    sphr_cntr_dstnc, cnstrnt_sphr_rds, ...
    dstnc_mtrcs{1}, dstnc_mtrcs{2});

problem_decomp = struct(...
    'x0', x0, ...
    'Aineq', actv_edg_orthgnl_cnstrnt.', ...
    'bineq', lnr_inqlty_rhs, ...
    'lb', [], ...
    'ub', [], ...
    'nonlcon', sphr_isos_cnstrnt , ...
    'objective', dstnc_mtrcs{1}, ...
    'solver', 'fmincon', ...
    'options', cndt_vrtx_optns);

% problem_decomp.objective = dstnc_objctvs_decomp{1};
% problem_decomp.nonlcon   = sphr_isos_cnstrnt;
% problem_decomp.x0        = x0;

%get the vertex coordinates that solve the constrained minimization problem
[new_vrtx_crdnts min_val exit_flag] = fmincon(problem_decomp);

%\/ DEBUG \/
%residual of constraint functions
%sphr_err = ...
%    (sphr_cntr_dstnc(new_vrtx_crdnts) - cnstrnt_sphr_rds)/cnstrnt_sphr_rds
%isos_err = ...
%    (dstnc_mtrcs{1}(new_vrtx_crdnts) - dstnc_mtrcs{2}(new_vrtx_crdnts))...
%    / (dstnc_mtrcs{1}(new_vrtx_crdnts) + dstnc_mtrcs{2}(new_vrtx_crdnts))
%/\ DEBUG /\

%dbg_h2 = plot3(tri_vert_coords(1), tri_vert_coords(2), tri_vert_coords(3), 'r*')
%delete(dbg_h1);
%delete(dbg_h2);
