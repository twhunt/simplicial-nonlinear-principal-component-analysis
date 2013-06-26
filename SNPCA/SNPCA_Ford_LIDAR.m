function [] = SNPCA_Ford_LIDAR()


lidar_data_file_path = 'Scan1000.mat';
load(lidar_data_file_path);

%crop scene
x_extnt = [-20 20];
y_extnt = [-20 20];
%x_extnt = [-10 15];
%y_extnt = [-10 15];


in_extnt = ...
    x_extnt(1) <= SCAN.XYZ(1, :) & SCAN.XYZ(1,:) <= x_extnt(2) ...
    & y_extnt(1) <= SCAN.XYZ(2, :) & SCAN.XYZ(2,:) <= y_extnt(2);

SCAN.XYZ = SCAN.XYZ(:, in_extnt);

%coordinates of initial point target
intl_pnt_crdnts = [9.46; -6.05; -.81]; %van
%intl_pnt_crdnts = [-9.13; .42; -.31]; %tree
%intl_pnt_crdnts = [-9.518; -.001183; -.2114];
%intl_pnt_crdnts = [11.34; -4.92; -1.21];
%intl_pnt_crdnts = [-2.56; 5.51; -2.45];
%intl_pnt_crdnts = [8.87; -6.48; -.66];
%intl_pnt_crdnts = [6.67; -1.18; .25];

dsplcmnt_crdnts = zeros(size(SCAN.XYZ));
for k=1:3
    dsplcmnt_crdnts(k, :) = SCAN.XYZ(k, :) - intl_pnt_crdnts(k);
end
dsplcmnt_lngth_sqrd = sum(dsplcmnt_crdnts.^2, 1);
[min_dsplcmnt_sqrd min_indx] = min(dsplcmnt_lngth_sqrd);
intl_pt_ind = min_indx;

disp(['Number surface data points: ' num2str(size(SCAN.XYZ,2))]);


SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .4;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = intl_pt_ind;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .5;
SNPCA_params.non_adj_tri_dist_tol          = .5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-6;
SNPCA_params.plot_frqncy                   = 1;
SNPCA_params.rtn_mtrx                      = eye(3);
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = false;
SNPCA_params.max_num_restarts              = 50;

[...
    vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals...
    ] ... 
    = SNPCA_interleaved_main(SCAN.XYZ, SNPCA_params);

