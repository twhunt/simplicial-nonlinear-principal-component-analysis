function [] = SNPCA_sphere()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate n random points that are distributed uniformly on the surface of
%the sphere
seed = 0;
rnd_strm = RandStream('mt19937ar','seed', seed);
RandStream.setGlobalStream(rnd_strm);
n   = 5e3; %number of data 

sphr_rds  = 4;

disp(['Sphere radius = ' num2str(sphr_rds)])
sphr_crdnts = gen_sphere_data_pts(sphr_rds, n);

%dimension of the original (before reduction) state space
state_space_dmnsn = 3;

[rtn_mtrx R] = qr(rand(state_space_dmnsn));
%rtn_mtrx = eye(3);
srfc_crdnts = rtn_mtrx(:,1:3)*sphr_crdnts;

%disp(['Initial search radius:' num2str(del)])
disp(['Number surface data points: ' num2str(n)])
%disp(['Noise magnitude: ' num2str(noise_magnitude)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
SNPCA_params.rtn_mtrx                      = rtn_mtrx;
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = true;
SNPCA_params.max_num_restarts              = 0;



[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);

