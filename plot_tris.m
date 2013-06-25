function tris_hndl = plot_tris(tris, tri_coords_x, tri_coords_y, tri_coords_z)

% colors =  mean(z(tris),2);
% colors = ((1+rd)/rd)*colors - rd;

tris_hndl = trisurf(tris,tri_coords_x,tri_coords_y,tri_coords_z);
%set alpha channel so triangulated surface is translucent
%alpha(.8)
alpha(1)
set(tris_hndl, 'FaceColor', [.4 .4 .4]);
%colormap(bone);
