function [tri_vrtx_inds, edg_vrtx_inds] = sew_seams_decomp4(...
    vrtx_crdnts, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    invbl_edg_inds, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    new_tri_max_edg_lngth, ...
    non_adj_tri_dist_tol, ...
    emprcl_drctn_crrltn_eval_bias, ...
    rtn_mtrx, ...
    plot_hndls, ...
    plot_frqncy)


vrtx_crdnts_3D = [];

[dmnsn num_vrtxs] = size(vrtx_crdnts);

dstnc_objctvs_decomp = cell(1,2);

num_cycls = 1;
while num_cycls > 0   
        
    tic
        
    frnt_cycl_edg_inds = frnt_cycl_basis(edg_vrtx_inds, tri_vrtx_inds);
    num_cycls          = numel(frnt_cycl_edg_inds);

    
    %remove all simple cycles of front edges that contain all inviable
    %front edges
    edg_is_invbl  = @(edg_ind) any(edg_ind == invbl_edg_inds);
    cycl_is_invbl = false(1, num_cycls);
    for cycl_i=1:num_cycls
       
        cycl_is_invbl(cycl_i) = ...
            all(arrayfun(edg_is_invbl, frnt_cycl_edg_inds{cycl_i}));
                
    end
    
    frnt_cycl_edg_inds(cycl_is_invbl) = [];
    num_cycls = numel(frnt_cycl_edg_inds);
    
    if num_cycls == 0
       
        disp('Exiting seam sewing stage');
        
    end
    
    edg_fifo_stack = [];
    for cycl_i=1:num_cycls
        
        edg_fifo_stack = [edg_fifo_stack; frnt_cycl_edg_inds{cycl_i}(:)];
        
    end
    
    srch_box_edg_lngths    = zeros(dmnsn, 1);
    srch_box_edg_lngths(:) = 1.5*new_tri_max_edg_lngth;
    num_to_find = 50;

    
    
    [tri_vrtx_inds, edg_vrtx_inds, invbl_edg_inds] ...
        = ...
        new_tris_frm_frnt_edg(...
        vrtx_crdnts, ...
        edg_fifo_stack, ...
        tri_vrtx_inds, ...
        edg_vrtx_inds, ...
        P_dmnt_egnvctrs, ...
        P_dmnt_egnvals, ...
        emprcl_drctn_crrltn_eval_bias, ...
        invbl_edg_inds, ...
        non_adj_tri_dist_tol, ...
        new_tri_max_edg_lngth, ...
        srch_box_edg_lngths,...
        num_to_find, ...
        vrtx_crdnts_3D);
    

    toc
    
    %save ruh_roh;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edg_edg_dot edg_edg_shrd_vrtx] = all_frnt_cycl_dots(...
    vrtx_crdnts, ...
    edg_vrtx_inds, ...
    frnt_cycls)

num_smpl_cycls = numel(frnt_cycls);

edg_edg_dot = cell(1, numel(frnt_cycls));
edg_edg_shrd_vrtx = cell(1, numel(frnt_cycls));

for cycl_i=1:num_smpl_cycls
    
    [edg_edg_dot{cycl_i} edg_edg_shrd_vrtx{cycl_i}] = ...
        frnt_cycl_dots(...
        vrtx_crdnts, ...
        edg_vrtx_inds, ...
        frnt_cycls{cycl_i});
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edg_edg_dot edg_edg_shrd_vrtx] = frnt_cycl_dots(...
    vrtx_crdnts, ...
    edg_vrtx_inds, ...
    frnt_cycl)


num_edgs = numel(frnt_cycl);
dmnsn    = size(vrtx_crdnts, 1);
 
edg_edg_dot       = zeros(1, num_edgs);
edg_edg_shrd_vrtx = zeros(1, num_edgs);



edg_sqnc_inds     = [0 0];
edg_vrtxs         = [0 0; 0 0];
non_shrd_vrtx_ind = [0 0];
unt_edg_vctr      = zeros(dmnsn, 2);

%recomputes unit vectors twice as many times as neccessary. kiss
for edg_sqnc_i=1:num_edgs
    
    edg_sqnc_inds(1) = edg_sqnc_i;
    if edg_sqnc_i ~= num_edgs
        
        edg_sqnc_inds(2) = edg_sqnc_i + 1;
        
    else
        
        edg_sqnc_inds(2) = 1;
        
    end
    
    for k=1:2

        edg_vrtxs(k, :) = edg_vrtx_inds(frnt_cycl(edg_sqnc_inds(k)), 1:2);
    
    end    
    
    tmp_shrd_vrtx_ind = intersect(edg_vrtxs(1,:), edg_vrtxs(2,:));
    if numel(tmp_shrd_vrtx_ind) ~= 1
       
        error('Number of shared vertices in adjacent front edges is not 1');
        
    end
    
    for k=1:2
       
        is_not_shrd = tmp_shrd_vrtx_ind ~= edg_vrtxs(k, :);
        tmp_not_shrd_vrtx_ind = edg_vrtxs(k, is_not_shrd);
        
        if numel(tmp_not_shrd_vrtx_ind) ~= 1
           
            error('Number of nonshared vertices in front edge is not 1')
            
        end
        
        non_shrd_vrtx_ind(k) = tmp_not_shrd_vrtx_ind;
        
        unt_edg_vctr(:, k) = ...
            vrtx_crdnts(:, non_shrd_vrtx_ind(k)) ...
            - vrtx_crdnts(:, tmp_shrd_vrtx_ind);
        
        unt_edg_vctr(:,k) = (1/norm(unt_edg_vctr(:,k)))*unt_edg_vctr(:,k);
        
    end
    
    edg_edg_dot(edg_sqnc_i) = unt_edg_vctr(:,1).'*unt_edg_vctr(:,2);
    edg_edg_shrd_vrtx(edg_sqnc_i) = tmp_shrd_vrtx_ind;
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nrst_vrtx_inds, sqrd_dstncs] = frnt_edg_nrst_vrtx_inds(...
    vrtx_crdnts, ...
    frnt_vrtx_inds, ...
    actv_edg_ind, ...
    edg_vrtx_inds, ...
    dstnc_objctvs_decomp, ...
    srch_box_edg_lngths, ...
    num_to_find)


[dmnsn num_vrtxs] = size(vrtx_crdnts);

actv_edg_vrtx_inds   = edg_vrtx_inds(actv_edg_ind, 1:2);
actv_edg_vrtx_crdnts = vrtx_crdnts(:, actv_edg_vrtx_inds);

frnt_vrtx_crdnts     = vrtx_crdnts(:, frnt_vrtx_inds);


%for each active edge vertex, find all vertices in the search box
extnts = zeros(dmnsn, 2);
hlf_srch_box_edg_lngths = .5*srch_box_edg_lngths(:);
is_in_actv_edg_vrtx_box = cell(1,2);
for vrtx_sqnc_i=1:2    
   
    extnts(:,1) = ...
        actv_edg_vrtx_crdnts(:,vrtx_sqnc_i) - hlf_srch_box_edg_lngths;
    extnts(:,2) = ...
        actv_edg_vrtx_crdnts(:,vrtx_sqnc_i) + hlf_srch_box_edg_lngths;

    is_in_actv_edg_vrtx_box{vrtx_sqnc_i} = ...
        in_box(frnt_vrtx_crdnts, extnts);

end

%logical indexing into frnt_edg_vrtcs
is_in_srch_rgn = ...
    is_in_actv_edg_vrtx_box{1} | is_in_actv_edg_vrtx_box{2};

%take out the actice edge vertices from consideration
is_in_srch_rgn(actv_edg_vrtx_inds(1) == frnt_vrtx_inds) = false;
is_in_srch_rgn(actv_edg_vrtx_inds(2) == frnt_vrtx_inds) = false;

num_in_srch_rgn = sum(is_in_srch_rgn);
in_srch_rgn_vrtx_inds = find(is_in_srch_rgn);
in_srch_rgn_vrtx_inds = frnt_vrtx_inds(in_srch_rgn_vrtx_inds);
actv_edgv_frnt_vrtx_dstnc_sqrd = zeros(num_in_srch_rgn, 2);


for vrtx_sqnc_i=1:num_in_srch_rgn
   
    vrtx_ind = in_srch_rgn_vrtx_inds(vrtx_sqnc_i);
    
    %     dbg_h = plot3(...
    %         vrtx_crdnts(1, vrtx_ind), ...
    %         vrtx_crdnts(2, vrtx_ind), ...
    %         vrtx_crdnts(3, vrtx_ind), 'r*')
    %
    %     delete(dbg_h);
    
    for k=1:2

        %         dlta_crdnts = vrtx_crdnts(:, vrtx_ind) - actv_edg_vrtx_crdnts(:,k);
        %         %actv_edgv_frnt_vrtx_dstnc_sqrd(vrtx_sqnc_i, k) = ...
        %         %    dlta_crdnts.'*(emprcl_lcl_drctn_crrltn_mtrx{k}*dlta_crdnts);
        %         actv_edgv_frnt_vrtx_dstnc_sqrd(vrtx_sqnc_i, k) = ...
        %             dstnc_objctvs_decomp{k}(dlta_crdnts);
        
        %dlta_crdnts = vrtx_crdnts(:, vrtx_ind) - actv_edg_vrtx_crdnts(:,k);
        %actv_edgv_frnt_vrtx_dstnc_sqrd(vrtx_sqnc_i, k) = ...
        %    dlta_crdnts.'*(emprcl_lcl_drctn_crrltn_mtrx{k}*dlta_crdnts);
        actv_edgv_frnt_vrtx_dstnc_sqrd(vrtx_sqnc_i, k) = ...
            dstnc_objctvs_decomp{k}(vrtx_crdnts(:, vrtx_ind));

        
    end
    
end

sum_actv_edgv_frnt_vrtx_dstnc_sqrd = ...
    sum(actv_edgv_frnt_vrtx_dstnc_sqrd, 2);

[sum_actv_edgv_frnt_vrtx_dstnc_sqrd sort_inds] = ...
    sort(sum_actv_edgv_frnt_vrtx_dstnc_sqrd, 'ascend');

num_to_retrn = min(num_to_find, numel(sum_actv_edgv_frnt_vrtx_dstnc_sqrd));

%nrst_vrtx_inds = zeros(num_to_retrn, 1);
%sqrd_dstncs    = zeros(num_to_retrn, 2);


nrst_vrtx_inds = in_srch_rgn_vrtx_inds(sort_inds(1:num_to_retrn));
sqrd_dstncs = actv_edgv_frnt_vrtx_dstnc_sqrd(sort_inds(1:num_to_retrn), :);

%nrst_vrtx_inds = nrst_vrtx_inds(end:-1:1);
%sqrd_dstncs    = sqrd_dstncs(end:-1:1,:);

end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tri_vrtx_inds, edg_vrtx_inds, invbl_edg_inds] = ...
    new_tris_frm_frnt_edg(...
    vrtx_crdnts, ...
    edg_fifo, ...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    emprcl_drctn_crrltn_eval_bias, ...
    invbl_edg_inds, ...
    non_adj_tri_dist_tol, ...
    new_tri_max_edg_lngth, ...
    srch_box_edg_lngths, num_to_find, ...
    vrtx_crdnts_3D)

%global vrtx_crdnts_3D;

%not really the max, but space for edg_fifo won't have to be reallocated if
%number of items in fifo buffer doesn't exceed MAX_EDGS
%INTL_FIFO_LNGTH = 64;

% INTL_NEW_TRIS_LGNTH = 32;
% INTL_NEW_EDGS_LNGTH = 32;
INTL_INVBLE_EDG_IND_LNGTH = 64;

% new_edg_vrtx_inds = edg_vrtx_inds;
% new_tri_vrtx_inds = tri_vrtx_inds;

%edg_fifo          = zeros(1, INTL_FIFO_LNGTH);
%edg_is_xstng      = false(1, INTL_FIFO_LNGTH);
% new_tri_vrtx_inds = zeros(INTL_NEW_TRIS_LGNTH, 3);
% new_edg_vrtx_inds = zeros(INTL_NEW_EDGS_LNGTH, 3);



ttl_num_new_edgs = 0;
ttl_num_new_tris = 0;

fifo_ind = numel(edg_fifo);
%edg_fifo(fifo_ind) = intl_frnt_edg_ind;

cndt_new_edgs = zeros(2,3);
is_exstng_edg = false(1,2);

edg_lngths = zeros(1,2);

tri_vrtx_inds = sort(tri_vrtx_inds, 2);

%frnt_cycl_edg_inds = frnt_cycl_basis(edg_vrtx_inds, tri_vrtx_inds);
dstnc_objctvs_decomp = cell(1,2);


is_frnt_edg = ...
    edge_blngs_to_xctly_one_tri(edg_vrtx_inds, tri_vrtx_inds);

%initialize is_frnt_edg to indicate that all edges are front edges
%entries in is_frnt_edg gets get set to correct values at top of while loop
%is_frnt_edg = true(size(edg_vrtx_inds,1), 1);

while fifo_ind ~= 0
    
    
    %dbg_new_tri  = zeros(1,3);
    dbg_new_edgs = zeros(2,3);

    crrnt_edg_is_frnt_edg = false;
    crrnt_edg_is_invbl    = true;
    stack_is_empty        = false;
    
        
    while (~crrnt_edg_is_frnt_edg || crrnt_edg_is_invbl ...
            || no_srfc_pnts_near_vrtxs) && ~stack_is_empty

        actv_edg_ind = edg_fifo(fifo_ind);
        edg_fifo(fifo_ind) = 0; %not strictly neccesary to zero out edge index
        %in array, but helpful for debugging
        fifo_ind = fifo_ind - 1

        stack_is_empty        = fifo_ind == 0;
        crrnt_edg_is_invbl    = any(actv_edg_ind == invbl_edg_inds);
        crrnt_edg_is_frnt_edg = is_frnt_edg(actv_edg_ind);

        no_srfc_pnts_near_vrtxs = ...
            isempty(P_dmnt_egnvctrs{edg_vrtx_inds(actv_edg_ind, 1)}) ...
            || isempty(P_dmnt_egnvctrs{edg_vrtx_inds(actv_edg_ind, 2)});
        
        if no_srfc_pnts_near_vrtxs
           
            invbl_edg_inds(end+1) = actv_edg_ind;
            
        end

        
    end
    
    if (~crrnt_edg_is_frnt_edg || crrnt_edg_is_invbl ...
            || no_srfc_pnts_near_vrtxs)
        %break out of outer while loop
       continue; 
    end    
       
    
    %unique list of all vertices belonging to a front edge
    frnt_vrtx_inds = unique(edg_vrtx_inds(is_frnt_edg,1:2));
    
    actv_edg_vrtx_inds = edg_vrtx_inds(actv_edg_ind, 1:2);
    actv_edg_vrtx_inds = actv_edg_vrtx_inds(:).';
    actv_edg_vrtx_inds = sort(actv_edg_vrtx_inds);
    
    actv_edg_crdnts = vrtx_crdnts(:, actv_edg_vrtx_inds);       
    
    for k=1:2        
        
        dstnc_objctvs_decomp{k} = @(x) dstnc_objctv_decomp(...
            x-actv_edg_crdnts(:,k), ...
            P_dmnt_egnvctrs{edg_vrtx_inds(actv_edg_ind, k)}, ...
            P_dmnt_egnvals{edg_vrtx_inds(actv_edg_ind, k)}, ...
            emprcl_drctn_crrltn_eval_bias);
        
    end

    
    %actv_edg_crdnts_3D = vrtx_crdnts_3D(:, actv_edg_vrtx_inds);
    %     dbg_actv_edg_hndl = plot3(...
    %        actv_edg_crdnts_3D(1, :), ...
    %        actv_edg_crdnts_3D(2, :), ...
    %        actv_edg_crdnts_3D(3, :), 'r-o');
    %
    %tmp_emprcl_lcl_drctn_crrltn_mtrx = ...
    %    {emprcl_lcl_drctn_crrltn_mtrx{actv_edg_vrtx_inds}};
    
    [nrst_vrtx_inds, sqrd_dstncs] = frnt_edg_nrst_vrtx_inds(...
        vrtx_crdnts, ...
        frnt_vrtx_inds, ...
        actv_edg_ind, ...
        edg_vrtx_inds, ...
        dstnc_objctvs_decomp, ...
        srch_box_edg_lngths, num_to_find);
    
        
    if isempty(nrst_vrtx_inds)
    
        invbl_edg_inds(end+1) = actv_edg_ind;
        %there was not a nearby front edge vertex
        continue
        %mark this edge as inviable
    end


    %\/ DEBUG \/
    %     nrstv_vrtx_hndl = plot3(...
    %         vrtx_crdnts_3D(1, nrst_vrtx_inds), ...
    %         vrtx_crdnts_3D(2, nrst_vrtx_inds), ...
    %         vrtx_crdnts_3D(3, nrst_vrtx_inds), ...
    %         'b+');
    %     %/\ DEBUG /\

    
    %num elts in nrst_vrtx_inds may be less than num_to_find
    for k=1:numel(nrst_vrtx_inds)
        
        nrst_vrtx_ind = nrst_vrtx_inds(k);
        
        %disp(['squared distance=' num2str(sum(sqrd_dstncs(k,:)))])
        
        %determine if the new triangle already exists
        tmp_new_tri_vrtx_inds = sort([actv_edg_vrtx_inds nrst_vrtx_ind]);
        
        new_tri = ...
            ~any(...
            tmp_new_tri_vrtx_inds(1) == tri_vrtx_inds(:, 1) ...
            & tmp_new_tri_vrtx_inds(2) == tri_vrtx_inds(:, 2) ...
            & tmp_new_tri_vrtx_inds(3) == tri_vrtx_inds(:, 3));
        
        if new_tri
            
            %determine if edges connecting active edge to new vertex are 
            %too long
            nrst_vrtx_crdnts = vrtx_crdnts(:, nrst_vrtx_ind);
            for kk=1:2
                                                
                edg_lngths(kk) = ...
                    norm(actv_edg_crdnts(:, kk) - nrst_vrtx_crdnts);
                
            end
            
            new_tri = new_tri && all(edg_lngths < new_tri_max_edg_lngth);
                                    
            if new_tri
                 
                frnt_edg_inds = ...
                    all_frnt_edg_inds(edg_vrtx_inds, tri_vrtx_inds);
                
                [is_error error_struct] = cand_vert_error_tmp(...
                    true, ...
                    nrst_vrtx_ind, ...
                    nrst_vrtx_crdnts, ...
                    actv_edg_ind, ...
                    frnt_edg_inds, ...
                    tri_vrtx_inds, ...
                    edg_vrtx_inds, ...
                    vrtx_crdnts, ...
                    non_adj_tri_dist_tol, ...
                    1);
                
                new_tri = ~is_error;
            
            end

        end
        
        if new_tri
            
            break
            
        end
        
    end

        
    if new_tri
                        
                
        cndt_new_edgs(1, :) = ...
            [sort([actv_edg_vrtx_inds(1) nrst_vrtx_ind]) ...
            actv_edg_vrtx_inds(2)];
        
        cndt_new_edgs(2, :) = ...
            [sort([actv_edg_vrtx_inds(2) nrst_vrtx_ind]) ...
            actv_edg_vrtx_inds(1)];
        
        new_tri_exstng_edg_ind = [0 0];
        for k=1:2
            
            %check if edge is in list of edges that were existing when this
            %function was called

            %vertices in edg_vrtx_inds may not be sorted
            vrtxs_mtch = ...
                cndt_new_edgs(k, 1) == edg_vrtx_inds(:, 1) ...
                & cndt_new_edgs(k, 2) == edg_vrtx_inds(:, 2);
            
            is_exstng_edg(k) = any(vrtxs_mtch);
            if ~is_exstng_edg(k)

                vrtxs_mtch = ...
                    cndt_new_edgs(k, 1) == edg_vrtx_inds(:, 2) ...
                    & cndt_new_edgs(k, 2) == edg_vrtx_inds(:, 1);

                is_exstng_edg(k) = any(vrtxs_mtch);
                
            end

            if is_exstng_edg(k)
                tmp_ind = find(vrtxs_mtch);
                assert(numel(tmp_ind) == 1);
                new_tri_exstng_edg_ind(k) = tmp_ind;
            end
            
            %             is_exstng_edg(k) = ...
            %                 any(...
            %                 cndt_new_edgs(k, 1) == edg_vrtx_inds(:, 1) ...
            %                 & cndt_new_edgs(k, 2) == edg_vrtx_inds(:, 2)) ...
            %                 || ...
            %                 any(...
            %                 cndt_new_edgs(k, 1) == edg_vrtx_inds(:, 2) ...
            %                 & cndt_new_edgs(k, 2) == edg_vrtx_inds(:, 1));
            
        end                
        
        num_new_edgs = 2 - sum(is_exstng_edg);
        
        ttl_num_new_tris = ttl_num_new_tris + 1;
        
        
        %         tri_vrtx_inds(end+1, :) = ...
        %             [sort(actv_edg_vrtx_inds) nrst_vrtx_ind];
        tri_vrtx_inds(end+1, :) = ...
            sort([actv_edg_vrtx_inds(:)' nrst_vrtx_ind]);
        
        dbg_new_tri = tri_vrtx_inds(end,:);

        
        %REMOVE this if statement \/
        %if any(~is_exstng_edg)
        %REMOVE this if statement /\
           
            for k=1:2
                
                if ~is_exstng_edg(k)
                    
                    % \/ DEBUG \/
                    %if size(edg_vrtx_inds,1) == 1061
                    %    warning('about to add a screwy edge')
                    %end
                    % /\ DEBUG /\
                    
                    edg_vrtx_inds(end+1, :) = cndt_new_edgs(k, :);
                    
                    fifo_ind = fifo_ind + 1;
                    %push the newly generated edge onto the stack
                    %its index is equal to the numer of edges in
                    %edg_vrtx_info
                    edg_fifo(fifo_ind) = size(edg_vrtx_inds, 1);
                    
                    dbg_new_edgs(k, :) = cndt_new_edgs(k, :);

                end
                
            end            
            
        %end
        ttl_num_new_edgs = ttl_num_new_edgs + num_new_edgs;
        
        %frnt_cycl_edg_inds = frnt_cycl_basis(edg_vrtx_inds, tri_vrtx_inds);

        is_frnt_edg(end:(end+num_new_edgs)) = true;
        
        %Update is front edge booleans
        %The new triangle always shares an edge with the active edge
        %Any existing edge that is part of the new triangle is no longer a
        %front edge
        is_frnt_edg(actv_edg_ind) = false;
        is_frnt_edg(new_tri_exstng_edg_ind(is_exstng_edg)) = false;
                
    else
        
        disp('Did not add a triangle')
        invbl_edg_inds(end+1) = actv_edg_ind;

    end
    
end

disp(['added ' num2str(ttl_num_new_edgs) ' edges and ' num2str(ttl_num_new_tris) ' triangles'])

if ttl_num_new_tris > 0
   warning('Marking all edges as valid. Likely inefficient!') 
   invbl_edg_inds = [];
end

end