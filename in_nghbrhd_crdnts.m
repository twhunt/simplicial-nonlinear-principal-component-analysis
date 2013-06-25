function is_in_ngbrhd = in_nghbrhd_crdnts(srfc_pt_crdnts, cntr_crdnts, rds)

dmnsn = size(srfc_pt_crdnts,1);

in_extnt = false(size(srfc_pt_crdnts));
for k=1:dmnsn
   
    in_extnt(k,:) = ...
        cntr_crdnts(k,:) - rds <= srfc_pt_crdnts(k,:) ...
        & srfc_pt_crdnts(k,:) <= cntr_crdnts(k,:) + rds;
    
end
    
is_in_ngbrhd   = all(in_extnt, 1);
num_pts_in_box = sum(is_in_ngbrhd);

crdnts_diff = zeros(dmnsn, num_pts_in_box);
for k=1:dmnsn
    
    crdnts_diff(k, :) = ...
        srfc_pt_crdnts(k, is_in_ngbrhd) - cntr_crdnts(k);
    
end

dstnc_sqrd       = sum(crdnts_diff.^2, 1);
in_box_not_in_ngbrhd = dstnc_sqrd > rds*rds;

tmp_inds = find(is_in_ngbrhd);
tmp_inds = tmp_inds(in_box_not_in_ngbrhd);

is_in_ngbrhd(tmp_inds) = false;




% dbg_crdnts = srfc_pt_crdnts(:, is_in_ngbrhd);
% dbg_crdnts_diff = zeros(dmnsn, size(dbg_crdnts,2));
% for k=1:dmnsn
%     
%     dbg_crdnts_diff(k, :) = ...
%         dbg_crdnts(k, :) - cntr_crdnts(k);
%     
% end
% 
% dbg_dstncs = sqrt(sum(dbg_crdnts_diff.^2, 1))