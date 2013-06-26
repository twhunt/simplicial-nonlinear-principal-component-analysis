function [] = SNPCA_sphere()

seed = 0; 
rnd_strm = RandStream('mt19937ar','seed', seed);
RandStream.setGlobalStream(rnd_strm);

state_space_dmnsn = 50;

SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = 1.5;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .3;
SNPCA_params.non_adj_tri_dist_tol          = .5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = .75*SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-4;
SNPCA_params.plot_frqncy                   = 1;
[SNPCA_params.rtn_mtrx, tmp_R]             = qr(rand(state_space_dmnsn));
%SNPCA_params.rtn_mtrx                      = eye(3);
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = true;
SNPCA_params.max_num_restarts              = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate n random points that are distributed uniformly on the surface of
%the sphere
n           = 5e3; %number of data 
sphr_rds    = 4;
sphr_crdnts = gen_sphere_data_pts(sphr_rds, n);
srfc_crdnts = SNPCA_params.rtn_mtrx(:,1:3)*sphr_crdnts;

%disp(['Initial search radius:' num2str(del)])
disp(['Number surface data points: ' num2str(n)])
disp(['Sphere radius = ' num2str(sphr_rds)])
%disp(['Noise magnitude: ' num2str(noise_magnitude)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);

