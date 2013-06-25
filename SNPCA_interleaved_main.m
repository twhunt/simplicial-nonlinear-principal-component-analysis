function [vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals]...
    = SNPCA_interleaved_main(srfc_crdnts, SNPCA_params)


MAX_NUM_VRTCS = 2048;

SAVE_FREQUENCY = 1;
%emprcl_drctn_crrltn_mtrx{vrtx_ind} holds the empirical local directin
%correlation matrix at the vertex associated with vrtx_ind
%emprcl_drctn_crrltn_mtrx = cell(1, MAX_NUM_VRTCS);

num_srfc_data_pnts = size(srfc_crdnts, 2);

%empirical covariance matrix P = P_factor*P_factor'
%P_factor        = cell(1, MAX_NUM_VRTCS);
%use map with two keys (one key string from catenated ieee to hex function)
P_dmnt_egnvctrs = cell(1, MAX_NUM_VRTCS);
P_dmnt_egnvals  = cell(1, MAX_NUM_VRTCS);


%store minimizer computed by constrained optimization problem so we only
%compute it once
%the coordinates of the candidate vertex generated from vertices with
%indices v1 and v2 are stored at intl_mnmzr_crdnts{v1,v2} and
%intl_mnmzr_crdnts{v2,v1}
%intl_mnmzr_crdnts = cell(MAX_NUM_VRTCS); %space inefficient!
intl_mnmzr_crdnts = {};

INTL_EDG_FIFO_LNGTH = 1024;
edg_fifo = zeros(1, INTL_EDG_FIFO_LNGTH);

srfc_pnt_is_vbl_intl = true(size(srfc_crdnts, 2), 1);

%surf_pt_blngs_to_tri:
%does the surface point belong to a triangle (boolean)
%if so, what triangle (triangles?) does it belong to?
%when placing a triangle vertex candidate, we may not want to consider
%surface points that near the surface of an existing triangle
% surf_pt_blngs_to_tri(num_surf_pts) = ...
%     struct('belongs', false(num_surf_pts,1), 'tri_inds', []);


%\/ Generate initial triangulation \/
%gen_init_tris2 generates an edge and two vertices so that the edge and the
%two vertices form the two initial triangles
[intl_edg_srfc_inds, intl_vrtx_inds] = gen_init_tris2(...
    srfc_crdnts, ...
    SNPCA_params.intl_pt_ind, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.cnstrnt_rad_fac*SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
    SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
    SNPCA_params.nnz_egnvals);

num_intl_tris = numel(intl_vrtx_inds);

if num_intl_tris == 0
   
    error('Error generating initial triangulation.')
    
else
    
    vrtx_ind_to_srfc_pt_ind = ...
        [intl_edg_srfc_inds(:).' intl_vrtx_inds(:).'];
    srfc_pt_ind_to_vrtx_ind = ...
        containers.Map(vrtx_ind_to_srfc_pt_ind, ...
        1:numel(vrtx_ind_to_srfc_pt_ind));
    
    vrtx_crdnts = srfc_crdnts(:, vrtx_ind_to_srfc_pt_ind);

    switch num_intl_tris
        case 2
        edg_vrtx_inds = [...
            1 2 3; ...
            2 3 1; ...
            1 3 2; ...
            2 4 1; ...
            1 4 2];
        
        tri_vrtx_inds = [...
            1 2 3; ...
            1 2 4];

        %push front edges of initial triangulation onto the edge stack
        edg_fifo(1:4) = [2 3 4 5];
        edg_fifo_ind  = 4;

        case 1
        
        edg_vrtx_inds = [...
            1 2 3; ...
            2 3 1; ...
            1 3 2];
        
        tri_vrtx_inds = [1 2 3];
        
        %push front edges of initial triangulation onto the edge stack
        edg_fifo(1:3) = [1 2 3];
        edg_fifo_ind  = 3;
        
            
    end
    
end
%/\ Generate initial triangulation /\


%prv_num_tris = size(tri_vrtx_inds, 1);

%surf_pt_is_tri_vert(ind) indicates whether the surface point indexed by
%ind is also a triangle vertex
%surf_pt_is_tri_vert is deprecated
%use srfc_pt_ind_to_vrtx_ind.isKey(srfc_pt_ind) where srfc_pt_ind is the
%index of the surface point into the vrtx_crds_x, etc.
%surf_pt_is_tri_vert = false(num_surf_pts,1);


%num_surf_pts = size(srfc_crdnts, 2);
data_dmnsn = size(srfc_crdnts, 1);

cndt_vrtx_info = new_cndt_vrtx_info(data_dmnsn);

%consider at most 10 candidate triangles
% cand_tri_vrtx_info(10) = cndt_vrtx_info;
% for k=1:9
%     cand_tri_vrtx_info(k) = cndt_vrtx_info;
% end


plot_hndls = new_plot_hndls();
if SNPCA_params.plot_frqncy > 0
    %plot surface points and initial triangulation if plotting is turned on
    figure(1);
    [az, el] = view();
    clf
    hold off
    
    
    %surface points
    srfc_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3).'*srfc_crdnts;
    plot_hndls.surf_pts_hndl = plot_surf_pts(...
        srfc_crdnts_3D(1, :), ...
        srfc_crdnts_3D(2, :), ...
        srfc_crdnts_3D(3, :));
    clear srfc_crdnts_3D ;
    
    %warning('Surface data points are invisible')
    %set(plot_hndls.surf_pts_hndl, 'Visible', 'off');
    
    disp('Clearing axis ticks')
    set(gca, ...
        'XTick', [], ...
        'YTick', [], ...
        'ZTick', [])
    
    %use open gl renderer for speed
    set(gcf, 'renderer', 'opengl')
    
    %clear dists extr_inds;
    
    view(az, el);
    axis equal
    axis vis3d;
    
    hold on
    
    %     plot_hndls = plot_tris_actv_edg(...
    %         plot_hndls, tri_vrtx_inds, edg_vrtx_inds, [], ...
    %         vrtx_crdnts, NLPCA_params);
    %
    vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3).'*vrtx_crdnts;
    plot_hndls.tris_hndl = plot_tris(tri_vrtx_inds, ...
        vrtx_crdnts_3D(1,:), vrtx_crdnts_3D(2,:), vrtx_crdnts_3D(3,:));
    
    axis equal
end


%sav_cnt = 0;

%dbg_num_tris      = size(tri_vrtx_inds,1);
%dbg_prvs_num_tris = dbg_num_tris;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEGIN load old program state and replot
%   %load; subplot(2,1,1); hold on; axis square; axis equal;
% load('saved_data/nt_2.mat'); %all_frnt_cycls_are_dead = false; seed is
% %load
%  plot_progress = true;
% %  load('tris_data_1700.mat'); plot_progress = true; clf; hold on; axis equal
% % seam_edg_inds = [];
% % front_is_subset_of_seams = false;
% surf_pts_hndl           = plot3([], [], [], '');
% actv_edge_hndl          = plot3([], [], [], '');
% actv_edge_intr_hndl     = plot3([], [], [], '');
% cand_vert_and_surf_hndl = plot3([], [], [], '');
% front_hndls             = plot3([], [], [], '');
% tris_hndl = plot_tris(tri_vrtx_inds, vrtx_crds_x, vrtx_crds_y, vrtx_crds_z);
%END load old program state and replot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% profile off
% profile on

%\/ skip by calling load \/
[...
    vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals...
    ]...
    = ...
    advancing_front_main_loop(...
    srfc_crdnts, ...
    edg_fifo, edg_fifo_ind, ...
    vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals, SNPCA_params, plot_hndls);

%mark all surface data points that coincide with vertices as inviable
%starting points
srfc_pnt_is_vbl_intl(vrtx_ind_to_srfc_pt_ind) = false;

if SAVE_FREQUENCY > 0
    disp('Saving initial advancing front run.');
    save([...
        'interleaved_data' ...
        filesep() ...
        'SNPCA_interleaved_run_' datestr(now, 30)]);
end


%save
%load
%load('interleaved_data/SNPCA_interleaved_run_20130623T151244.mat')

[tri_vrtx_inds, edg_vrtx_inds] = sew_seams_decomp4(...
    vrtx_crdnts, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    [], ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    SNPCA_params.new_tri_max_edg_lngth, ...
    SNPCA_params.non_adj_tri_dist_tol, ...
    SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
    SNPCA_params.rtn_mtrx, ...    
    plot_hndls, ...
    1);
%load
%save
%/\ skip by calling load /\

if SAVE_FREQUENCY > 0
    disp('Saving initial seam sewing run.');
    save([...
        'interleaved_data' ...
        filesep() ...
        'SNPCA_interleaved_run_' datestr(now, 30)]);
end

if SNPCA_params.plot_frqncy > 0
    vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtx_crdnts;
    plot_hndls.tris_hndl = plot_tris(...
        tri_vrtx_inds, ...
        vrtx_crdnts_3D(1,:), ...
        vrtx_crdnts_3D(2,:), ...
        vrtx_crdnts_3D(3,:));
end


fnd_vbl_intl_srfc_pnt = true;
plot_hndls_to_delete = {...
    'actv_edge_hndl', 'actv_edge_intr_hndl', 'front_hndls', ...
    'cand_vert_and_surf_hndl'};

restart_count = 0;

while fnd_vbl_intl_srfc_pnt ...
        && restart_count <= (SNPCA_params.max_num_restarts-1);

    
    restart_count = restart_count + 1;
    if mod(restart_count, SAVE_FREQUENCY) == 0
        disp(['Saving. Number of restarts: ' num2str(restart_count)]);
        save([...
            'interleaved_data' ...
            filesep() ...
            'SNPCA_interleaved_run_' datestr(now, 30)]);
    end
    

    [intl_tri_is_vbl, ...
    intl_edg_srfc_inds, ...
    intl_vrtx_inds,...
    srfc_pnt_is_vbl_intl] ...
    ...
    = intl_srfc_data_pnt(...
    ...
    srfc_crdnts, ...
    vrtx_crdnts, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    srfc_pnt_is_vbl_intl, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.cnstrnt_rad_fac*SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
    SNPCA_params.nnz_egnvals, ...
    SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
    SNPCA_params.non_adj_tri_dist_tol, ...
    num_srfc_data_pnts);

    fnd_vbl_intl_srfc_pnt = any(intl_tri_is_vbl);


    %restart if there's a data point that's far from the triangulation and
    %viable
    if fnd_vbl_intl_srfc_pnt
               
        % \/ add viable initial triangles to the existing triangulation \/
        max_vrtx_indx = max(tri_vrtx_inds(:));
        [...
            new_tri_vrtx_inds, new_edg_vrtx_inds, new_vrtx_inds, ...
            new_vrtx_crdnts, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind] ...
            = ...
            intl_cndt_tris_update(...
            intl_tri_is_vbl, intl_edg_srfc_inds, intl_vrtx_inds, ...
            max_vrtx_indx, srfc_crdnts, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind);

        tri_vrtx_inds = [tri_vrtx_inds; new_tri_vrtx_inds];
        edg_vrtx_inds = [edg_vrtx_inds; new_edg_vrtx_inds];
        vrtx_crdnts   = [vrtx_crdnts, new_vrtx_crdnts];
        
        
        num_new_edgs = size(new_edg_vrtx_inds, 1);
        num_edgs = size(edg_vrtx_inds, 1);
        edg_fifo(1:end) = 0;
        edg_fifo(1:num_new_edgs) = (num_edgs - num_new_edgs+1):num_edgs;
        edg_fifo_ind  = num_new_edgs;
        % /\ add viable initial triangles to the existing triangulation /\


        [...
            vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
            intl_mnmzr_crdnts, ...
            P_dmnt_egnvctrs, P_dmnt_egnvals...
            ]...
            = ...
            advancing_front_main_loop(...
            srfc_crdnts, ...
            edg_fifo, edg_fifo_ind, ...
            vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
            intl_mnmzr_crdnts, ...
            P_dmnt_egnvctrs, P_dmnt_egnvals, SNPCA_params, plot_hndls);
        
        %delete plot detritus
        if exist('plot_hndls', 'var') && isstruct(plot_hndls)
            for k=1:numel(plot_hndls_to_delete)
                
                if isfield(plot_hndls, plot_hndls_to_delete{k}) ...
                        && ~isempty(plot_hndls.(plot_hndls_to_delete{k})) ...
                        && ishandle(plot_hndls.(plot_hndls_to_delete{k}))
                    %plot_hndls.(plot_hndls_to_delete{k}) "dynamic field",
                    %access by string name of field
                    delete(plot_hndls.(plot_hndls_to_delete{k}));
                end
                
            end
            
        end

        
        [tri_vrtx_inds, edg_vrtx_inds] = sew_seams_decomp4(...
            vrtx_crdnts, ...
            tri_vrtx_inds, ...
            edg_vrtx_inds, ...
            [], ...
            P_dmnt_egnvctrs, ...
            P_dmnt_egnvals, ...
            SNPCA_params.new_tri_max_edg_lngth, ...
            SNPCA_params.non_adj_tri_dist_tol, ...
            SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
            SNPCA_params.rtn_mtrx, ...
            plot_hndls, ...
            1);

        %delete plot detritus
        if exist('plot_hndls', 'var') && isstruct(plot_hndls)
            for k=1:numel(plot_hndls_to_delete)
                
                if isfield(plot_hndls, plot_hndls_to_delete{k}) ...
                        && ~isempty(plot_hndls.(plot_hndls_to_delete{k})) ...
                        && ishandle(plot_hndls.(plot_hndls_to_delete{k}))
                    %plot_hndls.(plot_hndls_to_delete{k}) "dynamic field",
                    %access by string name of field
                    delete(plot_hndls.(plot_hndls_to_delete{k}));
                end
                
                
            end
            
        end
        
        if SNPCA_params.plot_frqncy > 0
            vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtx_crdnts;
            plot_hndls.tris_hndl = plot_tris(...
                tri_vrtx_inds, ...
                vrtx_crdnts_3D(1,:), ...
                vrtx_crdnts_3D(2,:), ...
                vrtx_crdnts_3D(3,:));
        end
        %delete_plot_lines(gca);
        
    end
    
end


%draw the tessellation and front one last time

if SNPCA_params.plot_frqncy > 0
    
    vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtx_crdnts;
    plot_hndls.tris_hndl = plot_tris(tri_vrtx_inds, ...
        vrtx_crdnts_3D(1,:), vrtx_crdnts_3D(2,:), vrtx_crdnts_3D(3,:));
    %set plot background color to white
    set(gcf, 'Color', [1 1 1])
    
end



% vrtx_crdnts_3D = NLPCA_params.rtn_mtrx(:,1:3).'*vrtx_crdnts;
%
% update_tri_surf(plot_hndls.tris_hndl, tri_vrtx_inds,...
%     vrtx_crdnts_3D(1,:), ...
%     vrtx_crdnts_3D(2,:), ...
%     vrtx_crdnts_3D(3,:));
%
% delete(plot_hndls.actv_edge_hndl(ishandle(plot_hndls.actv_edge_hndl)));
% delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
%     plot_hndls.actv_edge_intr_hndl)));
% delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
%     front_hndls)));
% delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
%     plot_hndls.cand_vert_and_surf_hndl)));

% front_hndls = ...
%     plot_front(frnt_cycl_edg_inds, frnt_cycl_ind, ...
%     edg_vrtx_inds, tri_vrtx_inds, front_info, ...
%     vrtx_crdnts_3D(1,:), ...
%     vrtx_crdnts_3D(2,:), ...
%     vrtx_crdnts_3D(3,:));


%save(['advancing_front_data' filesep() 'SNPCA_run_' datestr(now, 30)])

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%advancing_front_main_loop
function [...
    vrtx_crdnts, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    srfc_pt_ind_to_vrtx_ind, ...
    vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals]...
    = advancing_front_main_loop(...
    srfc_crdnts, ...
    edg_fifo, ...
    edg_fifo_ind, ...
    vrtx_crdnts, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    srfc_pt_ind_to_vrtx_ind, ...
    vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    SNPCA_params, plot_hndls)



%tell eigs that empirical covariance matrices are real symmetric
eigs_opts       = struct('issym', 1, 'isreal', 1);

%options for x0 unconstrained minimizer
%turn on user supplied gradient
x0_optns = optimset('fminunc');
x0_optns.GradObj = 'on';
%x0_optns.Display = 'iter';
%x0_optns.Diagnostics = 'on';
x0_optns.Display = 'off';


cndt_vrtx_optns = optimset('fmincon');
%num_surf_pts = size(srfc_crdnts, 2);
%data_dmnsn = size(srfc_crdnts, 1);

%distance as measured by induced metrics
dstnc_objctvs_decomp = cell(1,2);

%plotting variables
ordinate_index = 3; %1 for x, 2 for y, 3 for z. Used for triangle colo
local_color_map = bone();

while edg_fifo_ind > 0
    
    %sav_cnt = sav_cnt + 1;
    %save(['saved_data/nt_' num2str(sav_cnt) '.mat'])

    if SNPCA_params.plot_frqncy ~= 0
        delete(plot_hndls.actv_edge_hndl(...
            ishandle(plot_hndls.actv_edge_hndl)));
        delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
            plot_hndls.actv_edge_intr_hndl)));
        %delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
        %    front_hndls)));
        delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
            plot_hndls.cand_vert_and_surf_hndl)));
    end
    
    tic
    
    is_frnt_edg = ...
        edge_blngs_to_xctly_one_tri(edg_vrtx_inds, tri_vrtx_inds);

    %pop edges off the stack until the stack is empty, or popped edge is a
    %front edge
    crrnt_edg_is_frnt_edg = false;
    stack_is_empty = false;
    while ~stack_is_empty > 0 && ~crrnt_edg_is_frnt_edg
       
        actv_edge_ind          = edg_fifo(edg_fifo_ind);
        edg_fifo(edg_fifo_ind) = 0; %zero out edge index to aid debugging
        edg_fifo_ind           = edg_fifo_ind - 1;
        crrnt_edg_is_frnt_edg  = is_frnt_edg(actv_edge_ind);
        stack_is_empty         = edg_fifo_ind == 0;
        
    end
    
    if ~crrnt_edg_is_frnt_edg
        %condition at top of outer while loop will evaluate to false, and
        %advancing front stage will terminate
       continue; 
    end
        
    
    %edg_fifo(1:edg_fifo_ind)
    
    error_struct    = new_error_struct();
    cndt_vrtx_info  = new_cndt_vrtx_info(size(srfc_crdnts,1));
    
    vert_error      = false;
    
    cand_tri_ind    = NaN;
    
    %dbg_prvs_num_tris = dbg_num_tris;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: get indices and coordinates of active edge's vertices and
    % interior vertex
    actv_vert_inds        = edg_vrtx_inds(actv_edge_ind, 1:2);
    intr_vert_ind         = edg_vrtx_inds(actv_edge_ind, 3);
    actv_edge_vert_coords = vrtx_crdnts(:, actv_vert_inds(1:2));
    actv_intr_coords      = vrtx_crdnts(:, intr_vert_ind);
    actv_edg_vctr         = ...
        actv_edge_vert_coords(:,2) - actv_edge_vert_coords(:,1);
    actv_edg_lngth        = norm(actv_edg_vctr);
    actv_edg_mdpt_crdnts  = ...
        actv_edge_vert_coords(:,2) + actv_edge_vert_coords(:, 1);
    actv_edg_mdpt_crdnts  = .5*actv_edg_mdpt_crdnts;
    
    disp(['active edge length: ' num2str(actv_edg_lngth)]);
    % END: get indices and coordinates of active edge's vertices and
    % interior vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: plot tessellation, front, and active edge
    if mod(size(tri_vrtx_inds,1), SNPCA_params.plot_frqncy) == 0
        %delete(plot_hndls.actv_edge_hndl(...
        %    ishandle(plot_hndls.actv_edge_hndl)));
        %delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
        %    plot_hndls.actv_edge_intr_hndl)));
        %delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
        %    front_hndls)));
        %delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
        %    plot_hndls.cand_vert_and_surf_hndl)));
        
        vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtx_crdnts;
        
        plot_hndls = plot_tris_actv_edg(...
            plot_hndls, ...
            tri_vrtx_inds, edg_vrtx_inds, actv_edge_ind, ...
            vrtx_crdnts_3D);
        
        color_indices = tri_FaceVertexCData(tri_vrtx_inds,...
            vrtx_crdnts_3D(ordinate_index, :), ...
            size(local_color_map,1), ...
            min(vrtx_crdnts_3D(ordinate_index, :)),...
            max(vrtx_crdnts_3D(ordinate_index, :)));
        
        set(plot_hndls.tris_hndl, ...
            'CDataMapping', 'direct', ...
            'FaceVertexCData', color_indices(:));
        
    end
    % END: plot tessellation, front, and active edge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: place candidate vertex
    %try
    
    %DEBUG test
    %frnt_edg_adj_tris_driver(...
    %    edg_vrtx_inds, tri_vrtx_inds, frnt_cycl_edg_inds{frnt_cycl_ind}, actv_edge_ind, ...
    %    vrtx_crds_x, vrtx_crds_y, vrtx_crds_z)
    %end DEBUG test
            
    for vrtx_ind=edg_vrtx_inds(actv_edge_ind, 1:2)
        %                 k=find(isempty(emprcl_drctn_crrltn_mtrx(...
        %                     edg_vrtx_inds(actv_edge_ind, 1:2))))
        
        %calculate new direction empirical correlation matrix
        
        %get surface data points that are in the Euclidean
        %neighborhood of the current edge vertex
        is_in_nghbrhd = in_nghbrhd_crdnts(...
            srfc_crdnts, ...
            vrtx_crdnts(:, vrtx_ind), ...
            SNPCA_params.srch_rad_fac1*actv_edg_lngth);
        
        %take out surface data point that coincides with edge
        %vertex, so there isn't a zero vector when summing the
        %direction vector outer products
        is_in_nghbrhd(vrtx_ind_to_srfc_pt_ind(vrtx_ind)) = false;
        
        num_pts_in_ngbrhd = sum(is_in_nghbrhd);
        
        vert_error = num_pts_in_ngbrhd == 0;
        if vert_error
            %there are no surface data points in the search sphere
            %set vert_error to true
            warning('No surface data points in search sphere')
            break
            
        else
            
            %compute eigen decomposition of empirical direction
            %covariance matrices
            [P_dmnt_egnvctrs{vrtx_ind}, P_dmnt_egnvals{vrtx_ind}] = ...
                pnts_to_egn_dcmp(vrtx_crdnts(:, vrtx_ind), ...
                srfc_crdnts(:,is_in_nghbrhd), ...
                SNPCA_params.nnz_egnvals, ...
                SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
                eigs_opts);
            
        end
        
    end
        
            
    if ~vert_error
        %Set up metrics induced by perturbed inverse of empirical
        %covariance matrix
        
        %calculate constraint sphere radius as a weighted average of
        %the active edge length and the characteristic length
        cnstrnt_sphr_rds = ...
            (1-SNPCA_params.prfrd_cnstrnt_rds_wght)...
            *SNPCA_params.cnstrnt_rad_fac*actv_edg_lngth...
            + ...
            SNPCA_params.prfrd_cnstrnt_rds_wght...
            *SNPCA_params.prfrd_cnstrnt_rds;
        
        
        for k=1:2
            
            
            dstnc_objctvs_decomp{k} = @(x) dstnc_objctv_decomp2(...
                x-actv_edge_vert_coords(:,k), ...
                P_dmnt_egnvctrs{edg_vrtx_inds(actv_edge_ind, k)}, ...
                1./(P_dmnt_egnvals{edg_vrtx_inds(actv_edge_ind, k)}...
                + SNPCA_params.emprcl_drctn_crrltn_eval_bias), ...
                1/SNPCA_params.emprcl_drctn_crrltn_eval_bias);
            
            
        end
        
        dstnc_to_sphr_cntr = @(x) dstnc_objctv_decomp2(...
            x(:) - actv_edg_mdpt_crdnts, ...
            eye(numel(x)), ones(size(x(:),1),1), []);
        
        
        %dbg_x = ones(size(actv_edge_vert_coords(:,1)));
        %dbq_cmp = [dstnc_objctvs_decomp{1}(dbg_x) dstnc_objctvs_decomp_dbg{1}(dbg_x)]
        %/\ FROM gen_init_tris2 /\
        
        
        %Find initial point satisfying constraints
        [x0, eqlty_cnstrnt_val, x0_exit_flag] = cnstrnt_x0( ...
            dstnc_objctvs_decomp, ...
            dstnc_to_sphr_cntr, ...
            P_dmnt_egnvctrs(edg_vrtx_inds(actv_edge_ind, 1:2)), ...
            P_dmnt_egnvals(edg_vrtx_inds(actv_edge_ind, 1:2)), ...
            actv_edge_vert_coords, ...
            actv_intr_coords, ...
            actv_edg_mdpt_crdnts, ...
            cnstrnt_sphr_rds, ...
            SNPCA_params.emprcl_drctn_crrltn_eval_bias*ones(1,2), ...
            x0_optns);
        
        vert_error = x0_exit_flag <= 0;
        
        if vert_error
            warning(...
                ['Trouble finding initial point satisfying ' ...
                'equality constraints. fminunc exit flag: %d']...
                , x0_exit_flag);
        end
        
        
        if ~vert_error
            %find minimizer of constrained minimization problem
            [cand_tri_vert_coords, objctv_val, mnmzr_exit_flag] ...
                = gen_tri_vert_coords3_decomp(...
                actv_edge_vert_coords, ...
                actv_intr_coords, ...
                actv_edg_mdpt_crdnts, ...
                cnstrnt_sphr_rds, ...
                x0, ...
                dstnc_to_sphr_cntr, ...
                dstnc_objctvs_decomp, ...
                cndt_vrtx_optns);

                vert_error = mnmzr_exit_flag <= 0;

        end
        
        
        if vert_error
            warning(...
                'Problem finding minimizer. fmincon exit flag: %d', ...
                mnmzr_exit_flag);            
        end
        
        
        %save(...
        %    'pancake_subspace_info.mat', ...
        %    'srfc_crdnts', 'is_in_nghbrhd', ...
        %    'tmp_P_dmnt_egnvctrs', 'tmp_P_dmnt_egnvals', ...
        %    'actv_edge_vert_coords', 'actv_intr_coords', ...
        %    'NLPCA_params', ...
        %    'cand_tri_vert_coords');

        
    end
            
    %CLEAN_UP: don't use continue unless really necessary
    if vert_error
        continue
    end
        
    %candidate vertex was generated succesfuly
    %nudge the vertex to the nearest surface data point, and reject if
    %the distance between the original and nudged candidate vertices is
    %too great
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %START: nudge candidate vertex to nearest surface data point
    
    
    %search box edge length is twice maximum nudge distance
    %find data point whose mean distance in induced metrics is smallest
    %[cand_vert_nudge_dist nearest_surf_pt_ind] = ...
    %    nearest_non_tri_vert_surf_pt_Q(...
    %    cand_tri_vert_coords, ...
    %    srfc_crdnts, ...
    %    4*cand_vert_max_nudge_dist, ...
    %    dstnc_objctvs_decomp_cand_vert);
    
    if ~vert_error
        srch_rds = norm(cand_tri_vert_coords - actv_edg_mdpt_crdnts);
        
        is_in_nghbrhd = in_nghbrhd_crdnts(...
            srfc_crdnts, ...
            cand_tri_vert_coords, ...
            srch_rds);
        
        num_pts_in_ngbrhd = sum(is_in_nghbrhd);
        
        vert_error = num_pts_in_ngbrhd == 0;
    end
    
    if ~vert_error
        %\/ build metric based on surface data points near initial
        %candidate vertex \/
        [cand_vert_P_dmnt_egnvctrs, cand_vert_P_dmnt_egnvals] = ...
            pnts_to_egn_dcmp(...
            cand_tri_vert_coords, ...
            srfc_crdnts(:,is_in_nghbrhd), ...
            SNPCA_params.nnz_egnvals, ...
            SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
            eigs_opts);
        
        dstnc_objctv_decomp_cand_vert = @(x) dstnc_objctv_decomp(...
            x-cand_tri_vert_coords, ...
            cand_vert_P_dmnt_egnvctrs, ...
            cand_vert_P_dmnt_egnvals, ...
            SNPCA_params.emprcl_drctn_crrltn_eval_bias);
        %/\ build metric based on surface data points near initial
        %candidate vertex /\
        
        %find nearest vertex in Q metric associated with initial candidate
        %vertex
        [cand_vert_nudge_dist, nearest_surf_pt_ind] = ...
            nearest_srfc_pt_Q2(...
            cand_tri_vert_coords, ...
            srfc_crdnts, ...
            4*SNPCA_params.cand_vert_max_nudge_dist, ...
            dstnc_objctv_decomp_cand_vert);
                        
        cndt_vrtx_info.srfc_pt_ind = nearest_surf_pt_ind;
        
        
        %Was the minimizer nudged too far?
        vert_error = ...
            cand_vert_nudge_dist > SNPCA_params.cand_vert_max_nudge_dist;

        if vert_error
        
            warning(...
                ['Candidate vertex too far from surface data, ' ...
                'max distance "%e", actual distance "%e"'], ...
                SNPCA_params.cand_vert_max_nudge_dist, ...
                cand_vert_nudge_dist);
            
            %indicate candidate vertex placement error
            error_struct.err_id = 'NLPCA:too_much_nudge';
            
        
        end
        
        %Were there no surface data points near the minimizer?
        vert_error = vert_error || isempty(nearest_surf_pt_ind);
        
        if isempty(nearest_surf_pt_ind)
            warning([...
                'No surface points near minimizer in induced metric '...
                ' (this shouldnt happen, but is not fatal)'])
        end
        
    else
        
        %there are no surface data points in the search sphere
        %set vert_error to true
        warning('No surface data points in search sphere')
    
    end
    %END: nudge candidate vertex to nearest surface data point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: get tessellation info about candidate vertex        
    if ~vert_error
        frnt_edg_inds = find(is_frnt_edg);
        
        [cndt_vrtx_info.is_exstng_vrtx, ...
            cndt_vrtx_info.is_frnt_edg_vrtx,  ...
            cndt_vrtx_info.is_actv_edg_vrtx, ...
            cndt_vrtx_info.vrtx_ind] ...
            = srfc_pt_edg_info(...
            nearest_surf_pt_ind, actv_edge_ind, ...
            {frnt_edg_inds}, ...
            edg_vrtx_inds, srfc_pt_ind_to_vrtx_ind);
        
        cndt_vrtx_info.srfc_pt_ind      = nearest_surf_pt_ind;
        cndt_vrtx_info.crds             = ...
            srfc_crdnts(:, nearest_surf_pt_ind);
        
        
        %the active edge is not viable if the candidate vertex is a vertex
        %of the active edge (zero area), or it wasn't viable to begin with
        %\/ BUG? should be cndt_vrtx_info.is_frnt_edg_vrtx? \/
        %error_struct.actv_edg_is_vbl = ...
        %    ~is_actv_edg_vrtx && error_struct.actv_edg_is_vbl;
        %/\ BUG? /\
        %error_struct.actv_edg_is_vbl = ...
        %    ~cndt_vrtx_info.is_actv_edg_vrtx ...
        %&& error_struct.actv_edg_is_vbl;
        
        
        %vert_error = vert_error || ~error_struct.actv_edg_is_vbl;
        
        % END: get tessellation info about candidate vertex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % \/ Get list of non-shared vertices belonging to a front edge
        %adjacent to the active edge, and within the distance tolerance
        %specified by NLPCA_params.adj_vert_max_nudge_dist \/
        adjcnt_frnt_edg_cndt_vrtx_inds = cndt_adjcnt_vrtx_indxs(...
            actv_edge_ind, frnt_edg_inds, edg_vrtx_inds, ...
            cndt_vrtx_info.crds ,vrtx_crdnts, ...
            SNPCA_params.adj_vert_max_nudge_dist);
        % /\ Get list of non-shared vertices belonging to a front edge
        %adjacent to the active edge, and within the distance tolerance
        %specified by NLPCA_params.adj_vert_max_nudge_dist /\
        
        %candidate triangles are defined by the surface point nearest the
        %minimizer of the constrained optimization problem, and the
        %non-shared adjacent front edge vertices in
        %adjcnt_frnt_edg_cndt_vrtx_inds
        num_adj_cndt_vrtxs = numel(adjcnt_frnt_edg_cndt_vrtx_inds);
        num_cand_tris      = 1 + num_adj_cndt_vrtxs;
        
        
        %build list of info structs for candidate vertices
        %call cand_vert_error sequentially on entries of cand_tri_vrtx_info
        %so the nonadjacent candidate vertex is always checked last
        %add the adjacent existing candidate vertices in arbitrary order
        %(prefer smaller angles between adjacent active edges?)
        cand_tri_vrtx_info(num_adj_cndt_vrtxs+1) = cndt_vrtx_info;
        for k=1:num_adj_cndt_vrtxs
            
            cand_tri_vrtx_info(k) = cndt_vrtx_info;
            cand_tri_vrtx_info(k).is_exstng_vrtx   = true;
            cand_tri_vrtx_info(k).is_frnt_edg_vrtx = true;
            cand_tri_vrtx_info(k).is_actv_edg_vrtx = false;
            cand_tri_vrtx_info(k).vrtx_ind         = ...
                adjcnt_frnt_edg_cndt_vrtx_inds(k);
            cand_tri_vrtx_info(k).crds             = ...
                vrtx_crdnts(:, adjcnt_frnt_edg_cndt_vrtx_inds(k));
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %START: plot candidate vertex and nearest surface data point
        if mod(size(tri_vrtx_inds,1), SNPCA_params.plot_frqncy) == 0
            
            cand_tri_vert_coords_3D = ...
                SNPCA_params.rtn_mtrx(:,1:3).'*cand_tri_vert_coords;
            cndt_vrtx_info_crdnts_3D = ...
                SNPCA_params.rtn_mtrx(:,1:3).'*cndt_vrtx_info.crds(:);
            
            disp('Supressing plot of minimizer')
            plot_hndls.cand_vert_and_surf_hndl = plot3(...
                [cand_tri_vert_coords_3D(1) cndt_vrtx_info_crdnts_3D(1)], ...
                [cand_tri_vert_coords_3D(2) cndt_vrtx_info_crdnts_3D(2)], ...
                [cand_tri_vert_coords_3D(3) cndt_vrtx_info_crdnts_3D(3)], ...
                'x-r');
            
            drawnow
            
            %\/ animation \/
            %set(gcf, 'renderer', 'painter'); drawnow;
            %print('-depsc', ['advancing_front_animation/cand_vert_' num2str(sav_cnt)])
            if false
                disp('Generating advancing front animation')
                if ~exist('advancing_front_plot_count', 'var')
                    advancing_front_plot_count = 0;
                else
                    advancing_front_plot_count = ...
                        advancing_front_plot_count + 1;
                end
                print(...
                    '-dpng', '-r150', ...
                    ['animation/advancing_front/' ...
                    num2str(advancing_front_plot_count)]);
            end
            %/\ animation /\
        end
        %END: plot candidate vertex and nearest surface data point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % END: place candidate vertex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %if error_struct.actv_edg_is_vbl
        
        %candidate vertex is possibly viable, test if the candidate
        %triangle conflicts with an existing triangle
        
        nrby_frnt_edg_inds = nrby_frnt_edgs(...
            actv_edg_mdpt_crdnts, ...
            edg_vrtx_inds, ...
            tri_vrtx_inds, vrtx_crdnts, 2.5);
        
        cand_tri_cnflcts = false(1, num_cand_tris);
        for k=1:num_cand_tris
            %break out as soon as we detect a nonconflicting triangle.
            %the candidate triangeles should be sorted so the most
            %preferred triangle is at the
            
            [cand_tri_cnflcts(k), error_struct] = cand_vert_error_tmp(...
                cand_tri_vrtx_info(k).is_exstng_vrtx, ...
                cand_tri_vrtx_info(k).vrtx_ind, ...
                cand_tri_vrtx_info(k).crds, ...
                actv_edge_ind, ...
                nrby_frnt_edg_inds, ...
                tri_vrtx_inds, edg_vrtx_inds, ...
                vrtx_crdnts, ...
                SNPCA_params.non_adj_tri_dist_tol, 1);
            
            
            if ~cand_tri_cnflcts(k)
                %candidate triangle is acceptable, so don't test any
                %more.
                %alternative: pick best acceptable candidate triangle
                %instead of picking first acceptable candidate triangle
                cand_tri_ind = k;
                break
            end
            
        end
        
        %set vert_error to true if no candidate triangles, or all candidate
        %triangles conflicted
        vert_error = num_cand_tris == 0 || all(cand_tri_cnflcts);        
        % END: is the candidate vertex close enough to an existing
        % vertex? If so, then merge the two
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%
        
        
        %end
        
        %end
        
        
        %     dbg_num_tris  = size(tri_vrtx_inds, 1);
        %     added_tri = ~vert_error && error_struct.actv_edg_is_vbl;
        %
        %     assert(...
        %         (added_tri && dbg_num_tris ~= dbg_prvs_num_tris) ...
        %         || (~added_tri && dbg_num_tris == dbg_prvs_num_tris) );
        
        
    end
    
    if ~vert_error
        %Candidate triangle is acceptable. 
        
        if ~cand_tri_vrtx_info(cand_tri_ind).is_frnt_edg_vrtx
            %the candidate vertex is not an existing vertex and is
            %not a vertex of the front
            %one new vertex
            %two new edges
            
            num_new_edges = 2;
            
            new_vrtx_ind = size(vrtx_crdnts,2) + 1;
            
            vrtx_crdnts(:, new_vrtx_ind) = ....
                ...%cand_tri_vrtx_info(3).crds(:);
                cand_tri_vrtx_info(cand_tri_ind).crds(:);
            
            %associate the surface data point with the new vertex
            srfc_pt_ind_to_vrtx_ind(nearest_surf_pt_ind) = ...
                new_vrtx_ind;
            
            vrtx_ind_to_srfc_pt_ind(new_vrtx_ind) = nearest_surf_pt_ind;
            
            %the new triangle
            new_tri_vrtx_inds(1:2) = ...
                actv_vert_inds(1:2);
            new_tri_vrtx_inds(3) = new_vrtx_ind;
            
            %the 2 new edges
            new_edg_vrtx_inds(1,:) = ...
                [actv_vert_inds(1) ...
                new_vrtx_ind  actv_vert_inds(2)];
            
            new_edg_vrtx_inds(2,:) = ...
                [actv_vert_inds(2) ...
                new_vrtx_ind  actv_vert_inds(1)];
            
        else
            %the candidate vertex is an existing vertex in the
            %front
            
            %new_edgs = tri_vrtx_inds_to_edg_vrtx_inds(...
            %    cand_tri_vrtx_info(3).vrtx_ind, actv_edge_ind, ...
            %    edg_vrtx_inds, tri_vrtx_inds);
            
            new_edgs = tri_vrtx_inds_to_edg_vrtx_inds(...
                cand_tri_vrtx_info(cand_tri_ind).vrtx_ind, actv_edge_ind, ...
                edg_vrtx_inds, tri_vrtx_inds);
            
            
            
            num_new_edges = numel(new_edgs);
            
            for k=1:num_new_edges
                
                new_edg_vrtx_inds(k, :) = new_edgs{k};
                
            end
            
            new_tri_vrtx_inds(1:2) = actv_vert_inds(1:2);
            %new_tri_vrtx_inds(3)   = ...
            %    cand_tri_vrtx_info(3).vrtx_ind;
            new_tri_vrtx_inds(3)   = ...
                cand_tri_vrtx_info(cand_tri_ind).vrtx_ind;
            
        end
        
        tri_vrtx_inds(end+1,:) = new_tri_vrtx_inds;
        
        if num_new_edges == 2
            
            edg_vrtx_inds((end+1:end+2), :) = new_edg_vrtx_inds;
            
            %rjctd_cand_tri((end+1:end+2)) = false;
            
            for k=1:2
                
                edg_fifo_ind = edg_fifo_ind + 1;
                edg_fifo(edg_fifo_ind) = size(edg_vrtx_inds, 1)-(k-1);
                
            end
            
            
        elseif num_new_edges == 1
            
            edg_vrtx_inds(end+1, :) = new_edg_vrtx_inds(1, :);
            
            %rjctd_cand_tri(end+1) = false;
            
            edg_fifo_ind = edg_fifo_ind + 1;
            edg_fifo(edg_fifo_ind) = size(edg_vrtx_inds, 1);
            
            
        end
        
                        
    %else
                
        %num_new_edges                = 0;
        %error_struct.actv_edg_is_vbl = false;
        
    end
        
    toc
    
    %save ruh_roh_main
end




end %advancing_front_main_loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dstncs = pnt_vrtx_dstncs(...
    pnt_crdnts, ...
    vrtx_inds, ...
    vrtx_crdnts)

num_vrtxs = numel(vrtx_inds);
dstncs    = zeros(num_vrtxs,1);

if isempty(dstncs)
    return
end

dim_pnt_crdnts  = size(pnt_crdnts, 1);
dim_vrtx_crdnts = size(vrtx_crdnts, 1);

assert (dim_pnt_crdnts == dim_vrtx_crdnts);

dstncs = ( vrtx_crdnts(1, vrtx_inds) - pnt_crdnts(1)).^2;
for k=2:dim_vrtx_crdnts
    
    dstncs = dstncs + ( vrtx_crdnts(k, vrtx_inds) - pnt_crdnts(k) ).^2;
    
end

dstncs = sqrt(dstncs);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a front edge is near pnt_crdnts if a vertex of the front edge is within
%dstnc_tol of pnt_crdnts
function nrby_frnt_edg_inds = nrby_frnt_edgs(...
    pnt_crdnts, edg_vrtx_inds, tri_vrtx_inds, vrtx_crdnts, dstnc_tol)


%build list of vertices that belong to a front edge
frnt_edg_inds = all_frnt_edg_inds(edg_vrtx_inds, tri_vrtx_inds);
frnt_edg_vrtx_inds = edg_vrtx_inds(frnt_edg_inds, 1:2);

unq_frnt_edg_vrtx_inds = unique(frnt_edg_vrtx_inds(:));

%compute distance from the point with coordinates pnt_crdnts to front edge
%vertices
dstncs = pnt_vrtx_dstncs(pnt_crdnts, unq_frnt_edg_vrtx_inds, vrtx_crdnts);

dstnc_lt_tlrnc = dstncs < dstnc_tol;

nrby_frnt_vrtx_inds = unq_frnt_edg_vrtx_inds(dstnc_lt_tlrnc);

num_nrby_frnt_vrtxs = numel(nrby_frnt_vrtx_inds);

if num_nrby_frnt_vrtxs == 0
    
    nrby_frnt_edg_inds = [];
    
end

is_nrby_frnt_edg = ...
    nrby_frnt_vrtx_inds(1) == frnt_edg_vrtx_inds(:, 1) ...
    | nrby_frnt_vrtx_inds(1) == frnt_edg_vrtx_inds(:, 2);


for k=2:num_nrby_frnt_vrtxs
    
    is_nrby_frnt_edg = ...
        is_nrby_frnt_edg ...
        | ...
        nrby_frnt_vrtx_inds(k) == frnt_edg_vrtx_inds(:, 1) ...
        | nrby_frnt_vrtx_inds(k) == frnt_edg_vrtx_inds(:, 2);
    
end

nrby_frnt_edg_inds = frnt_edg_inds(is_nrby_frnt_edg);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



