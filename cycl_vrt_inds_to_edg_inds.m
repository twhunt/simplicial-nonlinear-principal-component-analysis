function cycl_edg_inds = cycl_vrt_inds_to_edg_inds(...
    cycl_vrt_inds, edg_vrt_inds)

num_cycl_vrt_inds = numel(cycl_vrt_inds);

if num_cycl_vrt_inds < 3
    error('Need at least 3 entries in array of vertex indices in cycle');
end

cycl_edg_inds = zeros(1, num_cycl_vrt_inds);

for k=1:num_cycl_vrt_inds-1
    
    crrnt_edg_vrt_ind1 = cycl_vrt_inds(k); 
    crrnt_edg_vrt_ind2 = cycl_vrt_inds(k+1);
    
    %the number of logical tests can be reduced if necessary
    %(find crrnt_edg_vrt_ind1 in edg_vrt_inds(:,1), then look at
    %corresponding locations in edg_vrt_inds(:,1) for crrnt_edg_vrt_ind2)
    is_edg = ...
        crrnt_edg_vrt_ind1 == edg_vrt_inds(:,1) ...
        & crrnt_edg_vrt_ind2 == edg_vrt_inds(:,2) ...
        | ...
        crrnt_edg_vrt_ind1 == edg_vrt_inds(:,2) ...
        & crrnt_edg_vrt_ind2 == edg_vrt_inds(:,1);
        
    crrnt_edg_ind = find(is_edg);

    if numel(crrnt_edg_ind) ~= 1
        error('Repeated edge in list of edge vertex indices')
    end
    
    cycl_edg_inds(k) = crrnt_edg_ind;
end

crrnt_edg_vrt_ind1 = cycl_vrt_inds(1);
%crrnt_edg_vrt_inds2 already set to cycl_vrt_inds(num_cycl_vrt_inds)

is_edg = ...
    crrnt_edg_vrt_ind1 == edg_vrt_inds(:,1) ...
    & crrnt_edg_vrt_ind2 == edg_vrt_inds(:,2) ...
    | ...
    crrnt_edg_vrt_ind1 == edg_vrt_inds(:,2) ...
    & crrnt_edg_vrt_ind2 == edg_vrt_inds(:,1);

crrnt_edg_ind = find(is_edg);

if numel(crrnt_edg_ind) ~= 1
    error('Repeated edge in list of edge vertex indices')
end

cycl_edg_inds(num_cycl_vrt_inds) = crrnt_edg_ind;