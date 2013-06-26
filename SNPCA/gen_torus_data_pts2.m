function trs_crdnts = gen_torus_data_pts2(...
    major_radius, minor_radius, num_pts)

uvs = [];
area = 4*pi^2*minor_radius*major_radius;

while size(uvs,2) <= num_pts

    unif_rand_grid = rand(3, num_pts);
    unif_rand_grid(1:2,:) = 2*pi*unif_rand_grid(1:2,:);
    unif_rand_grid(3, :) = ...
        (minor_radius*(minor_radius + major_radius))*unif_rand_grid(3,:);
    
    accpt = ...
        unif_rand_grid(3,:) ...
        <= ...
        minor_radius ...
        *(major_radius + minor_radius*cos(unif_rand_grid(2,:)));
        
    tmp_uvs = unif_rand_grid(:,accpt);
    uvs = [uvs, tmp_uvs];
end    

uvs = uvs(:, 1:num_pts);

trs_crdnts = zeros(3, num_pts);

trs_crdnts(1,:) = ...
    (major_radius + minor_radius*cos(uvs(2,:))).*cos(uvs(1,:));
trs_crdnts(2,:) = ...
    (major_radius + minor_radius*cos(uvs(2,:))).*sin(uvs(1,:));
trs_crdnts(3,:) = minor_radius*sin(uvs(2,:));

% plot(uvs(1,:), uvs(2,:), '.', 'MarkerSize', 1); axis equal;
% 
% plot3(trs_crdnts(1,:), trs_crdnts(2,:), trs_crdnts(3,:), '.', 'MarkerSize', 1);
% axis equal
% end