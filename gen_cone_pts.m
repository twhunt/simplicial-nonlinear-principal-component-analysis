function srfc_crdnts = gen_cone_pts(cone_radius, cone_height, num_pnts)

cone_surface_area = pi*cone_radius*sqrt(cone_height^2 + cone_radius^2);
disk_surface_area = pi*cone_radius^2;

ttl_srfc_area = cone_surface_area + disk_surface_area;

cone_prcntg = cone_surface_area/ttl_srfc_area;
disk_prcntg = disk_surface_area/ttl_srfc_area;

num_cone_pnts = floor(cone_prcntg*num_pnts);
num_disk_pnts = floor(disk_prcntg*num_pnts);

num_pnts_rmng = num_pnts - (num_cone_pnts + num_disk_pnts);
if num_pnts_rmng ~= 0
    
    num_cone_pnts = num_cone_pnts + 1;

end


cone_pdf_max = 1/pi;

prmtr_pnts = zeros(2,0);

%cone:
%prmtr_pnts(1, :) = radii (between 0 and 1)
%prmtr_pnts(1, :) = angles (between 0 and 2*pi)
while size(prmtr_pnts,2) <= num_cone_pnts

    unif_rand_grid = rand(3, num_cone_pnts);
    unif_rand_grid(2,:) = 2*pi*unif_rand_grid(2,:);
    unif_rand_grid(3, :) = cone_pdf_max*unif_rand_grid(3,:);
    
    accpt = ...
        unif_rand_grid(3,:) ...
        <= ...
        unif_rand_grid(1,:);
        
    %accpt = true(1,num_cone_pnts);
    
    tmp_prmtr_pnts = unif_rand_grid(:, accpt);
    prmtr_pnts = [prmtr_pnts, tmp_prmtr_pnts(1:2,:)];

end    

prmtr_pnts = prmtr_pnts(:, 1:num_cone_pnts);

srfc_crdnts = zeros(3, num_pnts);

srfc_crdnts(:, 1:num_cone_pnts) = [...
    cone_radius*prmtr_pnts(1,:).*cos(prmtr_pnts(2,:)); ...
    cone_radius*prmtr_pnts(1,:).*sin(prmtr_pnts(2,:)); ...
    cone_height*prmtr_pnts(1,:)];

% hold off
% plot3(...
%     srfc_crdnts(1, 1:num_cone_pnts), ...
%     srfc_crdnts(2, 1:num_cone_pnts), ...
%     srfc_crdnts(3, 1:num_cone_pnts), '.')
% axis equal


cone_pdf_max = 1/pi;
num_pnts_rmng = num_disk_pnts;

prmtr_pnts = zeros(2,0);

%cone:
%prmtr_pnts(1, :) = radii (between 0 and 1)
%prmtr_pnts(1, :) = angles (between 0 and 2*pi)
while size(prmtr_pnts,2) <= num_cone_pnts

    unif_rand_grid = rand(3, num_cone_pnts);
    unif_rand_grid(2,:) = 2*pi*unif_rand_grid(2,:);
    unif_rand_grid(3, :) = cone_pdf_max*unif_rand_grid(3,:);
    
    accpt = ...
        unif_rand_grid(3,:) ...
        <= ...
        unif_rand_grid(1,:);
        
    %accpt = true(1,num_cone_pnts);
    
    tmp_prmtr_pnts = unif_rand_grid(:, accpt);
    prmtr_pnts = [prmtr_pnts, tmp_prmtr_pnts(1:2,:)];

end    

prmtr_pnts = prmtr_pnts(:, 1:num_disk_pnts);

srfc_crdnts(1:2, (num_cone_pnts+1):end) = [...
    cone_radius*prmtr_pnts(1,:).*cos(prmtr_pnts(2,:)); ...
    cone_radius*prmtr_pnts(1,:).*sin(prmtr_pnts(2,:))];

srfc_crdnts(3, (num_cone_pnts+1):end) = cone_height;

% plot3(...
%     srfc_crdnts(1, (num_cone_pnts+1):end), ...
%     srfc_crdnts(2, (num_cone_pnts+1):end), ...
%     srfc_crdnts(3, (num_cone_pnts+1):end), '.')
% 
% 
% plot3(...
%     srfc_crdnts(1, :), ...
%     srfc_crdnts(2, :), ...
%     srfc_crdnts(3, :), '.')
% 
% axis equal