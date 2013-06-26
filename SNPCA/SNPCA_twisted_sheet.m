function [] = SNPCA_twisted_sheet()

seed = 1;
rnd_strm = RandStream('mt19937ar','seed', seed);
RandStream.setGlobalStream(rnd_strm);
num_pnts   = 1e4; %number of data 

%crdnts2D = randn(2, num_pnts);
midpnt_indx = round(num_pnts/2);
num_pnts_post_midpnt = num_pnts - midpnt_indx;

crdnts2D = zeros(2, num_pnts);
crdnts2D(:, 1:midpnt_indx) = rand(2, midpnt_indx);
crdnts2D(:, (midpnt_indx+1):end) = rand(2, num_pnts_post_midpnt);
crdnts2D(1, (midpnt_indx+1):end) = 1 + crdnts2D(1, (midpnt_indx+1):end);
crdnts2D(2, :) = crdnts2D(2, :) - .5;

max_twist_angls = zeros(4, 1);
max_twist_angls(:) = 2*pi;
max_twist_angls(end) = .5*max_twist_angls(end);

srfc_crdnts = twisted_plane2(crdnts2D, max_twist_angls);

% [coeff, score, latent] = pca(srfc_crdnts');
% latent

disp(['Number surface data points: ' num2str(num_pnts)])


SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .1;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .35;
SNPCA_params.non_adj_tri_dist_tol          = .5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = .75*SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-4;
SNPCA_params.plot_frqncy                   = 1;
SNPCA_params.rtn_mtrx                      = eye(size(srfc_crdnts, 1));
SNPCA_params.nnz_egnvals                   = numel(max_twist_angls)+2;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = false;
SNPCA_params.max_num_restarts              = 0;


[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);

