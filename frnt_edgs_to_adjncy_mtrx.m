function adjncy_mtrx = frnt_edgs_to_adjncy_mtrx(edg_vrtx_inds, frnt_edg_inds)

frnt_vrtx_inds = edg_vrtx_inds(frnt_edg_inds, 1:2);
frnt_vrtx_inds = sort(frnt_vrtx_inds, 2);

%\/ dense version \/
% adjncy_mtrx = false(max(frnt_vrtx_inds(:)));
% 
% tmp_lin_inds = ...
%     sub2ind(size(adjncy_mtrx), frnt_vrtx_inds(:,1), frnt_vrtx_inds(:,2));
% 
% adjncy_mtrx(tmp_lin_inds) = true;
% 
% tmp_lin_inds = ...
%     sub2ind(size(adjncy_mtrx), frnt_vrtx_inds(:,2), frnt_vrtx_inds(:,1));
% 
% adjncy_mtrx(tmp_lin_inds) = true;
%/\ dense version /\

adjncy_mtrx = sparse(...
    [frnt_vrtx_inds(:,1); frnt_vrtx_inds(:,2)], ...
    [frnt_vrtx_inds(:,2); frnt_vrtx_inds(:,1)], ...
    true(numel(frnt_vrtx_inds),1));