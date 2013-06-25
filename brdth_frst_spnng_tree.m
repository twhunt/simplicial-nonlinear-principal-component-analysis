function [spnng_tree bck_edgs1 bck_edgs2] = brdth_frst_spnng_tree(...
    adjncy_mtrx, strt_vrtx, max_bfs_queue_lngth)


%ensure that adjacency matrix is logical
if ~islogical(adjncy_mtrx)
    error('Adjacency matrix bust be of class logical');
end

spnng_tree = ...
    sparse(...
    1,1,false, ...
    size(adjncy_mtrx,1), size(adjncy_mtrx,1), ...
    nnz(adjncy_mtrx));

bck_edgs1 = [];
bck_edgs2 = [];

%graph vertices aren't necessarily numbered 1,2,...,n, so use an
%associative array
vrtx_vstd = containers.Map('KeyType', 'double', 'ValueType', 'logical');
vrtx_vstd(strt_vrtx) = true;

%first in first out buffer (queue)
%enqueue an element: fifo_arr(enqueue_ind)
%dequeue an element: fifo_arr(dequeue_ind)
fifo_arr    = zeros(1, max_bfs_queue_lngth);
fifo_arr(1) = strt_vrtx;

dequeue_ind = 1;
enqueue_ind = 2;

while dequeue_ind ~= enqueue_ind
    
    %dequeue the current vertex
    crrnt_vrtx = fifo_arr(dequeue_ind);
    
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
    cnnctd_vrtcs     = find(adjncy_mtrx(crrnt_vrtx, :));
    num_cnnctd_vrtcs = numel(cnnctd_vrtcs);
    
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
            
        elseif ~spnng_tree(crrnt_vrtx, cnnctd_vrtcs(k))
            %it's a back edge if it's not an edge in the spanning tree
            
            %edge connecting crrnt_vrtx and cnnctd_vrtcs(k) is a back edge
            %add back edge
            %edge is a back edge if it connects two vertices that have
            %already been visited
            %(the current vertex has been visited)
            if isempty(bck_edgs1)
                
                %bck_edgs = [crrnt_vrtx cnnctd_vrtcs(k)];
                bck_edgs1 = crrnt_vrtx;
                bck_edgs2 = cnnctd_vrtcs(k);
                
            elseif ~any( crrnt_vrtx == bck_edgs1 ...
                    & cnnctd_vrtcs(k) == bck_edgs2 ) ...
                    && ...
                    ~any( crrnt_vrtx == bck_edgs2 ...
                    & cnnctd_vrtcs(k) == bck_edgs1 )
                
                %bck_edgs = [bck_edgs; [crrnt_vrtx cnnctd_vrtcs(k)]];
                
                bck_edgs1(end+1) = crrnt_vrtx;
                bck_edgs2(end+1) = cnnctd_vrtcs(k);
            end
            
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUG
%bck_edgs
%spnng_tree
%DEBUG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
