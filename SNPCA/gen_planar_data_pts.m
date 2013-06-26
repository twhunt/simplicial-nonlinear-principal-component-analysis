function [crds_x crds_y] = gen_planar_data_pts(width, height, num_pts)

%generate random coordinates where 0 <= x <= 1, 0 <= y <= 1
crds_x = rand(num_pts,1);
crds_y = rand(num_pts,1);

%translate so mean is the origin
crd_mean_x = mean(crds_x);
crd_mean_y = mean(crds_y);

crds_x = width*(crds_x - crd_mean_x);
crds_y = height*(crds_y - crd_mean_y);

%plot(crds_x, crds_y, '.')