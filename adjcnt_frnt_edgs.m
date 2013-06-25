function adjcnt_frnt_edg_inds = adjcnt_frnt_edgs(...
    actv_edg_ind, frnt_edg_inds, edg_vrtx_inds)

actv_edg_vrtx_inds = edg_vrtx_inds(actv_edg_ind, :);

adjcnt_frnt_edg_inds = {[], []};

for k=1:2
    
    actv_edg_vrtx_ind = actv_edg_vrtx_inds(k);
    
    incdnt_edg_inds = vrtx_ind_to_incdnt_edg_inds(...
        actv_edg_vrtx_ind, frnt_edg_inds, edg_vrtx_inds);
    
    incdnt_frnt_edg_inds = intersect(frnt_edg_inds, incdnt_edg_inds);

    incdnt_frnt_edg_vrtx_inds = edg_vrtx_inds(incdnt_frnt_edg_inds, 1:2);
        
    %remove adjacent front edge belonging to active triangle (if it
    %exists)
    blngs_to_actv_tri = ...
        actv_edg_vrtx_inds(1) == incdnt_frnt_edg_vrtx_inds(:,1) ...
        & actv_edg_vrtx_inds(3) == incdnt_frnt_edg_vrtx_inds(:,2) ...
        | ...
        actv_edg_vrtx_inds(3) == incdnt_frnt_edg_vrtx_inds(:,1) ...
        & actv_edg_vrtx_inds(1) == incdnt_frnt_edg_vrtx_inds(:,2) ...
        | ...
        actv_edg_vrtx_inds(2) == incdnt_frnt_edg_vrtx_inds(:,1) ...
        & actv_edg_vrtx_inds(3) == incdnt_frnt_edg_vrtx_inds(:,2) ...
        | ...
        actv_edg_vrtx_inds(3) == incdnt_frnt_edg_vrtx_inds(:,1) ...
        & actv_edg_vrtx_inds(2) == incdnt_frnt_edg_vrtx_inds(:,2) ...
        | ...
        actv_edg_ind == incdnt_frnt_edg_inds;
        

    adjcnt_frnt_edg_inds{k} = incdnt_frnt_edg_inds(~blngs_to_actv_tri);
end

end