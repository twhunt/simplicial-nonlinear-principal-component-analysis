function [] = SNPCA_cone()


seed = 1; 
rnd_strm = RandStream('mt19937ar','seed', seed);
RandStream.setGlobalStream(rnd_strm);

state_space_dmnsn = 50;


SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .25;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .4;
SNPCA_params.non_adj_tri_dist_tol          = .5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-3;
SNPCA_params.plot_frqncy                   = 1;
[SNPCA_params.rtn_mtrx, tmp_R]             = qr(rand(state_space_dmnsn));
SNPCA_params.rtn_mtrx                      = eye(3);
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = false;
SNPCA_params.max_num_restarts              = 0;


n   = 1e4;  %number of data 

cone_height = 2; %horizontal torus radius
cone_radius = 1; %vertical torus radius
cone_crdnts = gen_cone_pts(cone_radius, cone_height, n);

disk_srfc_area = pi*cone_radius^2;
cone_srfc_area = pi*cone_radius*sqrt(cone_radius^2 + cone_height^2);
ttl_srfc_area  = cone_srfc_area + disk_srfc_area;
pnt_dnsty      = n/ttl_srfc_area;


srfc_crdnts = SNPCA_params.rtn_mtrx(:,1:3)*cone_crdnts;

disp(['Cone height: ' num2str(cone_height)])
disp(['Cone radius: ' num2str(cone_radius)])
disp(['state space dimension: ' num2str(state_space_dmnsn)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['Number surface data points: ' num2str(n)])
disp(['Surface point density: ' num2str(pnt_dnsty) ' points / area '])



[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);




