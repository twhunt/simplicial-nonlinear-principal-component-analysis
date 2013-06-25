function [] = SNPCA_creased_sheet()


state_space_dmnsn = 50;

seed = 0;
rnd_strm = RandStream('mt19937ar','seed', seed);
RandStream.setGlobalStream(rnd_strm);

SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .75;
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 1.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .5;
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
%Generate folded rectangular sheet
noise_magnitude = 0;
%noise_magnitude = .05;
%noise_magnitude = .01;

width = 4.25;
height = 11;
%the fold angle
fld_angl = .8;

n   = 1e4; %number of data points
[surf_coords_x surf_coords_y] = ...
    gen_planar_data_pts(width, height, round(.5*n));

%the left half of the sheet
surf_coords_x = surf_coords_x - .5*width;
surf_coords_z = zeros(size(surf_coords_x));

[plnr_crds_x plnr_crds_y] = ...
        gen_planar_data_pts(width, height, round(.5*n));

plnr_crds_x = plnr_crds_x + .5*width;
plnr_crds_z = zeros(size(plnr_crds_y));


for k=1:numel(plnr_crds_x)

    x_norm = abs(plnr_crds_x(k));
    plnr_crds_x(k) = x_norm*cos(fld_angl);
    plnr_crds_z(k) = x_norm*sin(fld_angl);

end

surf_coords_x = [surf_coords_x ; plnr_crds_x];
surf_coords_y = [surf_coords_y ; plnr_crds_y];
surf_coords_z = [surf_coords_z ; plnr_crds_z];

surf_coords_x = surf_coords_x + noise_magnitude*randn(size(surf_coords_x));
surf_coords_y = surf_coords_y + noise_magnitude*randn(size(surf_coords_y));
surf_coords_z = surf_coords_z + noise_magnitude*randn(size(surf_coords_z));


srfc_crdnts = SNPCA_params.rtn_mtrx(:,1:3)...
    *[surf_coords_x(:).'; surf_coords_y(:).'; surf_coords_z(:).'];

disp(['Number surface data points: ' num2str(n)])
disp(['Noise magnitude: ' num2str(noise_magnitude)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);

