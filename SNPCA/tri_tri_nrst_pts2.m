%calls quadprog instead of fmincon for 6x speed up
function [nrst_crds1, nrst_crds2, min_dstnc, exit_flag] ...
    = tri_tri_nrst_pts2(...
    tri1_crds1, tri1_crds2, tri1_crds3, ...
    tri2_crds1, tri2_crds2, tri2_crds3)


%tri1_edg1 = tri1_crds2(:) - tri1_crds1(:);
%tri1_edg2 = tri1_crds3(:) - tri1_crds1(:);

%tri2_edg1 = tri2_crds2(:) - tri2_crds1(:);
%tri2_edg2 = tri2_crds3(:) - tri2_crds1(:);

dT1 = [tri1_crds2(:) - tri1_crds1(:), tri1_crds3(:) - tri1_crds1(:)];
dT2 = [tri2_crds2(:) - tri2_crds1(:), tri2_crds3(:) - tri2_crds1(:)];

dT1dT1 = dT1.'*dT1;
dT2dT2 = dT2.'*dT2;
dT1_trnsps_dT2 = dT1.'*dT2;

d_crds = tri1_crds1 - tri2_crds1;

d_crds_lnght_sqrd = sum(d_crds.^2);

%points on the two triangles are encoded as:
%a1*tri1_edg1 + a2*tri1_edg2
%b1*tri2_edg1 + b2*tri2_edg2
%where:
%a1, a2, b1, b2 >= 0
%a1 + a2 <= 1
%b1 + b2 <= 1
%the variable edge_coefs is arranged as [a1; a2; b1; b2]


%linear constraints: A*edge_coefs <= b
A = ...
    [1  1  0  0; ...
    -1  0  0  0; ...
     0 -1  0  0; ...
     0  0  1  1; ...
     0  0 -1  0; ...
     0  0  0 -1];
 
b = [1; 0; 0; 1; 0 ;0];

dstnc_sqrd_hssn = 2*[dT1dT1 -dT1_trnsps_dT2; -dT1_trnsps_dT2.' dT2dT2];

dstnc_sqrd_grdnt = 2*[dT1'*d_crds; -dT2'*d_crds];

persistent options;
if isempty(options)
    options = optimset('quadprog');
    options = optimset(options, ...
    'Algorithm', 'interior-point-convex', ...
    'Display', 'off');

end



[mnmzg_coefs, min_dstnc, exit_flag] = quadprog(...
    dstnc_sqrd_hssn, dstnc_sqrd_grdnt, A, b, ...
    [], [], [], [], [], ...
    options);

if exit_flag < 0
    warning('fmincon problem in tri_tri_nrst_pts')
    disp(['exit_flag = ' num2str(exit_flag)]);
end

%objective function in quadprog has zero constant term, so add original
%constant term to the minimized distance between the two triangles
%we minimized the squared distance, so take a square root
min_dstnc = d_crds_lnght_sqrd + min_dstnc;
min_dstnc = sqrt(abs(min_dstnc));

% mnmzg_coefs(1)*tri1_edg1 + mnmzg_coefs(2)*tri1_edg2
% mnmzg_coefs(3)*tri2_edg1 + mnmzg_coefs(4)*tri2_edg2

nrst_crds1 = tri1_crds1(:) + dT1*mnmzg_coefs(1:2);

nrst_crds2 = tri2_crds1(:) + dT2*mnmzg_coefs(3:4);

return

%objective function is the squared distance between the two points
%associated with
% dstnc_sqrd = @(edge_coefs) ...
%     dot(...
%     orgn_diff ...
%     + edge_coefs(1)*tri1_edg1 + edge_coefs(2)*tri1_edg2 ...
%     - edge_coefs(3)*tri2_edg1 - edge_coefs(4)*tri2_edg2, ...
%     orgn_diff ...
%     + edge_coefs(1)*tri1_edg1 + edge_coefs(2)*tri1_edg2 ...
%     - edge_coefs(3)*tri2_edg1 - edge_coefs(4)*tri2_edg2);

% dstnc_sqrd = @(x) objctv(...
%     x, orgn_diff, dT1, dT2, dT1dT1, dT2dT2, dT1_trnsps_dT2);

dstnc_sqrd = @(x) objctv2(...
    x, tri1_crds1(:), dT1, tri2_crds1(:), dT2);


dstnc_sqrd_hssn = @(x, lambda) ...
    2*[dT1dT1 -dT1_trnsps_dT2; -dT1_trnsps_dT2.' dT2dT2];

%initial coefficients for minimization
intl_coefs = [.25; .25; .25; .25];  

options = optimset('fmincon');
options = optimset(options, ...
    'Algorithm', 'interior-point', ...
    'Display', 'notify-detailed', ...
    'GradObj', 'on', ...
    'Hessian', 'user-supplied', ...
    'HessFcn', dstnc_sqrd_hssn);

% options = optimset(options, ...
%     'Algorithm', 'interior-point', ...
%     'Display', 'notify-detailed', ...
%     'GradObj', 'on');

%[mnmzg_coefs min_dstnc] = fmincon(dstnc_sqrd, intl_coefs, A, b);
[mnmzg_coefs min_dstnc exit_flag] = ...
    fmincon(dstnc_sqrd,intl_coefs,A,b,[],[],[],[],[],options);

if exit_flag < 0
    warning('fmincon problem in tri_tri_nrst_pts')
    disp(['exit_flag = ' num2str(exit_flag)]);
end

%we minimized the squared distance, so take a square root
min_dstnc = sqrt(min_dstnc);

% mnmzg_coefs(1)*tri1_edg1 + mnmzg_coefs(2)*tri1_edg2
% mnmzg_coefs(3)*tri2_edg1 + mnmzg_coefs(4)*tri2_edg2

nrst_crds1 = tri1_crds1(:) + dT1*mnmzg_coefs(1:2);

nrst_crds2 = tri2_crds1(:) + dT2*mnmzg_coefs(3:4);




% dbg_tri1_crdnts = [tri1_crds1(:) tri1_crds2(:) tri1_crds3(:) tri1_crds1(:)];
% dbg_tri2_crdnts = [tri2_crds1(:) tri2_crds2(:) tri2_crds3(:) tri2_crds1(:)];
% 
% dbg_tri_hndl = plot3(...
%     dbg_tri1_crdnts(1, :), ...
%     dbg_tri1_crdnts(2, :), ...
%     dbg_tri1_crdnts(3, :), 'k-', ...
%     dbg_tri2_crdnts(1, :), ...
%     dbg_tri2_crdnts(2, :), ...
%     dbg_tri2_crdnts(3, :), 'b-');
% 
% dbg_mnmzr_hndl = plot3(...
%     nrst_crds1(1), ...
%     nrst_crds1(2), ...
%     nrst_crds1(3), '*k',...
%     nrst_crds2(1), ...
%     nrst_crds2(2), ...
%     nrst_crds2(3), '*b');
% 
% delete(dbg_tri_hndl);
% delete(dbg_mnmzr_hndl);


end