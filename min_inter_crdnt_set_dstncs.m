%[min_dstnc, min_indx] = min_inter_crdnt_set_dstncs(in_crdnts1, in_crdnts2)
%compute the shortest distance from each point in in_crdnts1 to a point in
%in_crdnts2
%min_dstnc(k) is the minimum distance from the point with coordinates
%in_crdnts(:, k) to a point in in_crdnts2
%min_indx(k) holds the index of the point in in_crdnts2 nearest the point
%with coordinates in_crdnts1(:, k)
function [min_dstnc, min_indx] = ...
    min_inter_crdnt_set_dstncs(in_crdnts1, in_crdnts2)

dmnsn = size(in_crdnts1, 1);
if dmnsn ~= size(in_crdnts2, 1)

    error('Dimension mismatch in coordinate sets')
    
end

num_pnts    = zeros(1,2);
num_pnts(1) = size(in_crdnts1, 2);
num_pnts(2) = size(in_crdnts2, 2);

min_dstnc = zeros(1, num_pnts(1));
min_indx  = zeros(1, num_pnts(1));

dsplcmnts  = zeros(dmnsn, num_pnts(2));
for pnt_k=1:num_pnts(1)

    for dmnsn_k=1:dmnsn

        dsplcmnts(dmnsn_k, :) = ...
            in_crdnts2(dmnsn_k, :) - in_crdnts1(dmnsn_k, pnt_k);
    
    end
        
    tmp_dstncs = sum(dsplcmnts.^2, 1);
    
    [min_dstnc(pnt_k), min_indx(pnt_k)] = min(tmp_dstncs);
    min_dstnc(pnt_k) = sqrt(min_dstnc(pnt_k));
    
end