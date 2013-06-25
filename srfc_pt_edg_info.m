function [is_vrtx, is_frnt_edg_vrtx, is_actv_edg_vrtx, vrtx_ind] = ...
    srfc_pt_edg_info(...
    srfc_pt_ind, actv_edg_ind, frnt_cycl_edg_inds, ...
    edg_vrtx_inds, srfc_pt_ind_to_vrtx_ind)


is_vrtx = srfc_pt_ind_to_vrtx_ind.isKey(srfc_pt_ind);

if is_vrtx

    %index to a vertex (triangle or edge)
    vrtx_ind = srfc_pt_ind_to_vrtx_ind(srfc_pt_ind);
    
    %is the vertex a vertex of the active edge
    is_actv_edg_vrtx = any(vrtx_ind == edg_vrtx_inds(actv_edg_ind, 1:2)); 
    
    if is_actv_edg_vrtx
        %the vertex is a vertex of the active edge, so it's a front edge
        %vertex
        is_frnt_edg_vrtx = true;
    else
        
        %is the vertex a vertex of a front edge?
        for k=1:numel(frnt_cycl_edg_inds)
            
            tmp_frnt_cycl_edg_inds = ...
                edg_vrtx_inds(frnt_cycl_edg_inds{k}, 1:2);
            
            is_frnt_edg_vrtx = ...
                any(vrtx_ind == tmp_frnt_cycl_edg_inds(1:end));
            
            assert(numel(is_frnt_edg_vrtx) == 1);
            
            if is_frnt_edg_vrtx
                break;
            end
            
        end
    end
        
    
else
    
    is_frnt_edg_vrtx = false;
    is_actv_edg_vrtx = false;
    vrtx_ind = NaN;
    
end

