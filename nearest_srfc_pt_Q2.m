function [min_Q_dstnc, ind_of_nearest] = nearest_srfc_pt_Q2(...
    pt_crdnts, ...
    surf_crdnts, ...
    srch_box_edg_lngth, ...
    dstnc_objctvs_decomp)


% find all surface coordinates inside a box centered at pt_crdnts, then
% compute distances in induced metric
half_srch_box_edg_lngth = .5*srch_box_edg_lngth;

in_extent = ...
    pt_crdnts(1) - half_srch_box_edg_lngth <= surf_crdnts(1,:) ...
    & ...
    surf_crdnts(1,:) <= pt_crdnts(1) + half_srch_box_edg_lngth;

for k=2:size(surf_crdnts,1)

    in_extent = in_extent & ...
        pt_crdnts(k) - half_srch_box_edg_lngth <= surf_crdnts(k,:) ...
        & ...
        surf_crdnts(k,:) <= pt_crdnts(k) + half_srch_box_edg_lngth;
        
end

cndt_inds = find(in_extent);
cndt_crdnts = surf_crdnts(:, in_extent);
Q_dstncs = zeros(1, numel(cndt_inds));

%first, compute squared distances in induced metric
for k=1:numel(cndt_inds)
   
    Q_dstncs(k) = dstnc_objctvs_decomp(cndt_crdnts(:,k));
        
    
end

%compute distance, instead of squared distance
%Q_dstncs = sqrt(Q_dstncs);

%Q_dstncs_mean = .5*sum(Q_dstncs, 1);

[min_Q_dstnc min_ind] = min(Q_dstncs);
min_Q_dstnc           = sqrt(min_Q_dstnc);
ind_of_nearest        = cndt_inds(min_ind);
