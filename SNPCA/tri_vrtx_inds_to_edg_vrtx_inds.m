function new_edgs = tri_vrtx_inds_to_edg_vrtx_inds(...
    cand_vrtx_ind, actv_edg_ind, ...
    edg_vrtx_inds, tri_vrtx_inds)

%build new edges when candidate vertex is an existing vertex
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(cand_vrtx_ind)
    disp('empty cand_vrtx_ind');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actv_edg_vrtx_inds = edg_vrtx_inds(actv_edg_ind, 1:2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START: generate list of all front edges
frnt_edg_inds = all_frnt_edg_inds(edg_vrtx_inds, tri_vrtx_inds);

actv_edg_incdnt_edg_inds  = cell(1, 2);
incdnt_edg_is_new_tri_edg = cell(1,2);
actv_edg_vrtx_is_intr     = false(1,2);

for k=1:2

    actv_edg_incdnt_edg_inds{k} = vrtx_ind_to_incdnt_edg_inds(...
        actv_edg_vrtx_inds(k), frnt_edg_inds, edg_vrtx_inds);
    
    actv_edg_incdnt_edg_inds{k} = setdiff(...
        actv_edg_incdnt_edg_inds{k}, actv_edg_ind);
    
    incdnt_edg_is_new_tri_edg{k} = ...
        false(1, numel(actv_edg_incdnt_edg_inds));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %START: determine if edge incident to the current active edge vertex
    %owns the candidate vertex
    
    for kk=1:numel(actv_edg_incdnt_edg_inds{k})
       
        incdnt_edg_is_new_tri_edg{k}(kk) = ...
            (cand_vrtx_ind ...
            == edg_vrtx_inds(actv_edg_incdnt_edg_inds{k}(kk), 1) ...
            && ...
            actv_edg_vrtx_inds(k) ...
            == edg_vrtx_inds(actv_edg_incdnt_edg_inds{k}(kk), 2) ) ...
            || ...
            (cand_vrtx_ind ...
            == edg_vrtx_inds(actv_edg_incdnt_edg_inds{k}(kk), 2) ...
            && ...
            actv_edg_vrtx_inds(k) ...
            == edg_vrtx_inds(actv_edg_incdnt_edg_inds{k}(kk), 1) );
        
    end

    assert(nnz(incdnt_edg_is_new_tri_edg{k}) <= 1);
    
    actv_edg_vrtx_is_intr(k) = any(incdnt_edg_is_new_tri_edg{k});
    %END: determine if edge incident to the current active edge vertex
    %owns the candidate vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

num_new_edgs = 2 - ...
    nnz(incdnt_edg_is_new_tri_edg{1}) - nnz(incdnt_edg_is_new_tri_edg{2});

switch num_new_edgs
    case 0
        %the new triangle is composed of 3 existing edges
        %this happens when the front cycle is composed of exactly three
        %front edges
        new_edgs = {};
        
    case 1
        
        %incdnt_adj_edg_ind = ...
        %    adj_frnt_edg_vrtx_inds(cand_vrtx_is_adj_edg_vrtx);
        
        %         is_cand_vrtx = ...
        %             cand_vrtx_ind ...
        %             == adj_frnt_edg_vrtx_inds(cand_vrtx_is_adj_edg_vrtx, 1:2);
        %
        %         intr_vrtx_ind = adj_frnt_edg_vrtx_inds(...
        %             cand_vrtx_is_adj_edg_vrtx, ~is_cand_vrtx);
        %
        %         edg_vrtx_ind = ...
        %             actv_edg_vrtx_inds(intr_vrtx_ind ~= actv_edg_vrtx_inds);
        
        new_edgs = {[...
            cand_vrtx_ind ...
            actv_edg_vrtx_inds(~actv_edg_vrtx_is_intr)...
            actv_edg_vrtx_inds(actv_edg_vrtx_is_intr) ]};
    case 2
        
        new_edgs = {...
            [cand_vrtx_ind actv_edg_vrtx_inds(1) actv_edg_vrtx_inds(2)],...
            [cand_vrtx_ind actv_edg_vrtx_inds(2) actv_edg_vrtx_inds(1)]};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START: get front cycle edges adjacent to the active edge
frnt_cycl_vrtx_inds = edg_vrtx_inds(frnt_cycl_edg_inds, 1:2);

is_adj_frnt_edg = false(size(frnt_cycl_vrtx_inds));

is_adj_frnt_edg(:, 1) = ...
    actv_edg_vrtx_inds(1) == frnt_cycl_vrtx_inds(:, 1) ...
    & actv_edg_vrtx_inds(2) ~= frnt_cycl_vrtx_inds(:, 2) ...
    | actv_edg_vrtx_inds(1) == frnt_cycl_vrtx_inds(:, 2) ...
    & actv_edg_vrtx_inds(2) ~= frnt_cycl_vrtx_inds(:, 1);

is_adj_frnt_edg(:, 2) = ...
    actv_edg_vrtx_inds(2) == frnt_cycl_vrtx_inds(:, 1) ...
    & actv_edg_vrtx_inds(1) ~= frnt_cycl_vrtx_inds(:, 2) ...
    | actv_edg_vrtx_inds(2) == frnt_cycl_vrtx_inds(:, 2) ...
    & actv_edg_vrtx_inds(1) ~= frnt_cycl_vrtx_inds(:, 1);

if nnz(is_adj_frnt_edg(:, 1)) ~= 1 || nnz(is_adj_frnt_edg(:, 2)) ~= 1
    error('Exactly one front edge must be incident to each vertex of the active edge')
end

adj_frnt_edg_inds = ...
    [frnt_cycl_edg_inds(is_adj_frnt_edg(:, 1)) ...
     frnt_cycl_edg_inds(is_adj_frnt_edg(:, 2))];
% END: get front cycle edges adjacent to the active edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
adj_frnt_edg_vrtx_inds = zeros(2, 2);
adj_frnt_edg_vrtx_inds(1, :) = edg_vrtx_inds(adj_frnt_edg_inds(1), 1:2);
adj_frnt_edg_vrtx_inds(2, :) = edg_vrtx_inds(adj_frnt_edg_inds(2), 1:2);

cand_vrtx_is_adj_edg_vrtx = ...
    [any(cand_vrtx_ind == adj_frnt_edg_vrtx_inds(1, :)) ...
     any(cand_vrtx_ind == adj_frnt_edg_vrtx_inds(2, :))];

num_new_edgs = 2 - nnz(cand_vrtx_is_adj_edg_vrtx);
 
switch num_new_edgs
    case 0
        %the new triangle is composed of 3 existing edges
        %this happens when the front cycle is composed of exactly three
        %front edges
        new_edgs = {};
        
    case 1
        
        %incdnt_adj_edg_ind = ...
        %    adj_frnt_edg_vrtx_inds(cand_vrtx_is_adj_edg_vrtx);
        
        is_cand_vrtx = ...
            cand_vrtx_ind ...
            == adj_frnt_edg_vrtx_inds(cand_vrtx_is_adj_edg_vrtx, 1:2);
        
        intr_vrtx_ind = adj_frnt_edg_vrtx_inds(...
            cand_vrtx_is_adj_edg_vrtx, ~is_cand_vrtx);
        
        edg_vrtx_ind = ...
            actv_edg_vrtx_inds(intr_vrtx_ind ~= actv_edg_vrtx_inds);
        
        new_edgs = {[cand_vrtx_ind edg_vrtx_ind intr_vrtx_ind]};
    case 2
        
        new_edgs = {...
            [cand_vrtx_ind actv_edg_vrtx_inds(1) actv_edg_vrtx_inds(2)],...
            [cand_vrtx_ind actv_edg_vrtx_inds(2) actv_edg_vrtx_inds(1)]};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end