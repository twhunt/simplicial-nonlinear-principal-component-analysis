function surf_pts_hndl = plot_surf_pts(surf_coords_x, surf_coords_y, surf_coords_z)

surf_pts_hndl = plot3(surf_coords_x,surf_coords_y,surf_coords_z,'.');
set(surf_pts_hndl, 'MarkerSize', 2);
set(surf_pts_hndl, 'MarkerFaceColor', [0 .8 .8]);
set(surf_pts_hndl, 'MarkerEdgeColor', [0 .8 .8]);