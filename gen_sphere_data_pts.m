function sphr_crdnts = gen_sphere_data_pts(radius, num_pts)


%num_accptd = 0;
thetas = [];
area = 4*pi*radius^2;
while size(thetas, 2) < num_pts

    unif_rand_grid = rand(3, num_pts);
    unif_rand_grid(1,:) = 2*pi*unif_rand_grid(1,:);
    unif_rand_grid(2,:) = pi*unif_rand_grid(2,:);
    unif_rand_grid(3,:) = radius^2*unif_rand_grid(3,:);
    
    accpt = unif_rand_grid(3,:) ...
        <= (radius^2*sin(unif_rand_grid(2,:)))/area;
    tmp_thetas = unif_rand_grid(1:2, accpt);
    thetas = [thetas, tmp_thetas];
end    

thetas = thetas(:, 1:num_pts);

cos_theta = cos(thetas(1,:));
sin_theta = sin(thetas(1,:));

cos_phi = cos(thetas(2,:));
sin_phi = sin(thetas(2,:));

sphr_crdnts = radius*[...
    cos_theta.*sin_phi; ...
    sin_theta.*sin_phi; ...
    cos_phi];


% plot(thetas(1,:), thetas(2,:), '.', 'MarkerSize', 1); axis equal;
% 
% plot3(sphr_crdnts(1,:), sphr_crdnts(2,:), sphr_crdnts(3,:), '.', 'MarkerSize', 1);
% axis equal
% end
% 

%plot3(tor_pts_x, tor_pts_y, tor_pts_z, '.')
%axis equal
