function color_indices = tri_FaceVertexCData(tri_vrtx_inds,...
                              ordinates, ...
                              color_map_max_index, range_min, range_max)
                          
%normalized_range will hold values between 0 and 1, which will then map to
%the discrete set 1, 2, ..., color_map_max_index
for k=size(tri_vrtx_inds,1):-1:1
    
    color_indices(k) = mean(ordinates(tri_vrtx_inds(k,:)));
    
end

color_indices = ...
    (color_map_max_index/(range_max - range_min))...
    *(color_indices - range_min);

 color_indices = ceil(color_indices);                                                   