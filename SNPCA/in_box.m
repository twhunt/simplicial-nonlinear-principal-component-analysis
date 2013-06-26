function is_in_box = in_box(clmn_crdnts, extnts)

in_extnt = false(size(clmn_crdnts));
for k=1:size(clmn_crdnts, 1)

    in_extnt(k, :) = extnts(k, 1) <= clmn_crdnts(k, :)...
        & clmn_crdnts(k, :) <= extnts(k, 2);

end

is_in_box = all(in_extnt, 1);
