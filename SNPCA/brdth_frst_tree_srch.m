function path = brdth_frst_tree_srch(...
    tree_adjncy_mtrx, strt_vrtx, end_vrtx, max_bfs_queue_lngth)


%ensure that adjacency matrix is logical
if ~islogical(tree_adjncy_mtrx)
    error('Adjacency matrix bust be of class logical');
end



spnng_tree = ...
    sparse(...
    1, 1, false, ...
    size(tree_adjncy_mtrx,1), size(tree_adjncy_mtrx,1), ...
    nnz(tree_adjncy_mtrx));


%graph vertices aren't necessarily numbered 1,2,...,n, so use an
%associative array
vrtx_vstd = containers.Map('KeyType', 'double', 'ValueType', 'logical');
vrtx_vstd(strt_vrtx) = true;

vrtx_tree_dpth = ...
    containers.Map('KeyType', 'double', 'ValueType', 'double');
vrtx_tree_dpth(strt_vrtx) = 1;

crrnt_tree_dpth = 1;

%first in first out buffer (queue)
%enqueue an element: fifo_arr(enqueue_ind)
%dequeue an element: fifo_arr(dequeue_ind)
fifo_arr    = zeros(1, max_bfs_queue_lngth);
fifo_arr(1) = strt_vrtx;

dequeue_ind = 1;
enqueue_ind = 2;

fnd_end_vrtx = false;

while dequeue_ind ~= enqueue_ind && ~fnd_end_vrtx
    
    
    %dequeue the current vertex
    crrnt_vrtx = fifo_arr(dequeue_ind);

    crrnt_tree_dpth = vrtx_tree_dpth(crrnt_vrtx) + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DEBUG
    %fifo_arr
    %fifo_arr(dequeue_ind) = 0;
    %DEBUG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %increment circular array index by 1
    if dequeue_ind < max_bfs_queue_lngth
        
        dequeue_ind = dequeue_ind + 1;
        
    else
        
        dequeue_ind = 1;
        
    end
    
    %list of all vertices connected to the current vertex
    cnnctd_vrtcs     = find(tree_adjncy_mtrx(crrnt_vrtx, :));
    num_cnnctd_vrtcs = numel(cnnctd_vrtcs);
    %cnnctd_vrtx_vstd = false;
    
    fnd_end_vrtx = any(end_vrtx == cnnctd_vrtcs);
    
    if fnd_end_vrtx
        
        spnng_tree(...
            sub2ind(size(spnng_tree), ...
            [crrnt_vrtx end_vrtx], ...
            [end_vrtx crrnt_vrtx])) ...
            = true;
        
        vrtx_tree_dpth(end_vrtx) = crrnt_tree_dpth;
    else
        for k=1:num_cnnctd_vrtcs
            
            %the connected vertex has been visited if it is in the
            %(associative) array cnnctd_vrtcs
            cnnctd_vrtx_vstd = vrtx_vstd.isKey(cnnctd_vrtcs(k));
            
            if ~cnnctd_vrtx_vstd
                
                %edge connecting current vertex and connected vertex has not
                %been traversed
                
                %enqueue all unvisited vertices connected to to the current
                %vertex
                %each of these vertices is the endpoing of an edge in the
                %spanning tree. the other edge vertex is the current vertex
                
                vrtx_vstd(cnnctd_vrtcs(k)) = true;
                fifo_arr(enqueue_ind) = cnnctd_vrtcs(k);
                
                vrtx_tree_dpth(cnnctd_vrtcs(k)) = crrnt_tree_dpth;
                
                spnng_tree(...
                    sub2ind(size(spnng_tree), ...
                    [crrnt_vrtx cnnctd_vrtcs(k)], ...
                    [cnnctd_vrtcs(k) crrnt_vrtx])) ...
                    = true;
                
                %increment circular array index by 1
                if enqueue_ind < max_bfs_queue_lngth
                    
                    enqueue_ind = enqueue_ind + 1;
                    
                else
                    
                    enqueue_ind = 1;
                    
                end
                
                if dequeue_ind == enqueue_ind
                    error('Queue filled up while building bfs spanning tree')
                end
                
            end
            
        end
    end
    
end


%generate the path connecting the start and end vertices, if it exists
if ~fnd_end_vrtx
    path = [];
else
    path = zeros(1, crrnt_tree_dpth);
    path(crrnt_tree_dpth) = end_vrtx;

    for k=crrnt_tree_dpth:-1:2

        cnnctd_vrtcs = find(spnng_tree(path(k), :));
        
        is_prnt = false;
        for kk=1:numel(cnnctd_vrtcs)
            is_prnt = ...
                vrtx_tree_dpth.isKey(cnnctd_vrtcs(kk)) ...
                && vrtx_tree_dpth(cnnctd_vrtcs(kk)) == k-1;
            if is_prnt 
                path(k-1) = cnnctd_vrtcs(kk);
                break;
            end
            
        end
    end
end
