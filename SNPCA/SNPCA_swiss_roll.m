function [] = SNPCA_swiss_roll()

seed = 0;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

state_space_dmnsn = 50;


SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .75;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .35;
SNPCA_params.non_adj_tri_dist_tol          = .5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-3;
SNPCA_params.plot_frqncy                   = 1;
[SNPCA_params.rtn_mtrx, tmp_R]                  = qr(rand(state_space_dmnsn));
%SNPCA_params.rtn_mtrx                      = eye(3);
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = false;
SNPCA_params.max_num_restarts              = 0;



num_pnts   = 1e4; %number of data 
min_angle = 0;
max_angle = 2*(2*pi);
crvtr_prmtr = .5;
width = 6;
swiss_roll_crdnts = ...
    swiss_roll(num_pnts, crvtr_prmtr, width, min_angle, max_angle);

srfc_crdnts = SNPCA_params.rtn_mtrx(:,1:3)*swiss_roll_crdnts;


disp(['Number surface data points: ' num2str(num_pnts)])
disp(['curvature factor:' num2str(crvtr_prmtr)])
disp(['min/max angle:' num2str(min_angle) ' ' num2str(max_angle)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);
