%plot_hndls = plot_tris_actv_edg(
%plot_hndls, tri_vrtx_inds, edg_vrtx_inds, actv_frnt_edg_ind, vrtx_crdnts,
%NLPCA_params
%pass in [] for actv_frnt_edg_ind to supress plotting of active edge and 
%interior vertex
function plot_hndls = plot_tris_actv_edg(...
    plot_hndls, ...
    tri_vrtx_inds, edg_vrtx_inds, actv_frnt_edg_ind, vrtx_crdnts)

%erase old plot data
%sometimes the graphics is invalid, and delete fails, so only
%delete graphics objects with valid handles
%it would be nice to know why the graphics handles are
%intermittently invalid
delete(plot_hndls.actv_edge_hndl(ishandle(plot_hndls.actv_edge_hndl)));
delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
    plot_hndls.actv_edge_intr_hndl)));
delete(plot_hndls.front_hndls(ishandle(...
    plot_hndls.front_hndls)));
delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
    plot_hndls.cand_vert_and_surf_hndl)));

%Dumb! pass in the 3D cooordinates!
%vrtx_crdnts_3D = NLPCA_params.rtn_mtrx(:,1:3)'*vrtx_crdnts;

%add new triangle to the existing surface
%without creating a whole new surface
update_tri_surf(plot_hndls.tris_hndl, tri_vrtx_inds,...
    vrtx_crdnts(1,:), vrtx_crdnts(2,:), vrtx_crdnts(3,:));

if ~isempty(actv_frnt_edg_ind)
    actv_tri_crdnts = ...
        vrtx_crdnts(:, edg_vrtx_inds(actv_frnt_edg_ind, :));
    %actv_tri_crdnts(:, 1:2) are the active edge coordinates
    %actv_tri_crdnts(:, 3) are the coordinates of the interior vertex of the
    %active edge
    
    %vertices of the active edge
    %disp('Supressing active edge vertices plot')
    plot_hndls.actv_edge_hndl = plot3(...
        actv_tri_crdnts(1, 1:2), ...
        actv_tri_crdnts(2, 1:2), ...
        actv_tri_crdnts(3, 1:2), 'ob');
    
    %interior vertex of the active edge
    %disp('Supressing active edge interior vertex plot')
    plot_hndls.actv_edge_intr_hndl = plot3(...
        actv_tri_crdnts(1, 3), ...
        actv_tri_crdnts(2, 3), ...
        actv_tri_crdnts(3, 3), 'sb');
end
%\/ animation \/
%set(gcf, 'renderer', 'painter'); drawnow;
%print('-depsc', ['advancing_front_animation/new_tri_' num2str(sav_cnt)])
%/\ animation /\

end