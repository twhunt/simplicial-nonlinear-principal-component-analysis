function is_in_ngbrhd = in_nghbrhd(srfc_pt_crdnts, cntr_ind, rds)

cntr_crdnts = srfc_pt_crdnts(:, cntr_ind);

is_in_ngbrhd = in_nghbrhd_crdnts(srfc_pt_crdnts, cntr_crdnts, rds);

return