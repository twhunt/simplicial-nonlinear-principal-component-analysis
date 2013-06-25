function crdnts = swiss_roll(num_pnts, crvtr, width, min_angl, max_angl)

%generate num_pts parameter pairs

angls_taus = [];
area = ...
    .5*crvtr...
    *(max_angl*sqrt(1+max_angl^2) - min_angl*sqrt(1-min_angl^2) ...
    + asinh(max_angl) - asinh(min_angl));

while size(angls_taus,2) <= num_pnts

    unif_rand_grid = rand(3, num_pnts);
    
    unif_rand_grid(1,:) = ...
        min_angl + (max_angl - min_angl)*unif_rand_grid(1,:);
    unif_rand_grid(2,:) = width*unif_rand_grid(2,:);
    unif_rand_grid(3,:) = sqrt(1+min_angl^2)*unif_rand_grid(3,:);
    
    accpt = unif_rand_grid(3,:) ...
        <= ...
        (crvtr*sqrt(1+unif_rand_grid(1,:).^2))/area;
        
    tmp_angls_taus = unif_rand_grid(:,accpt);
    angls_taus = [angls_taus, tmp_angls_taus];
end    

angls_taus = angls_taus(:, 1:num_pnts);

crdnts = zeros(3, num_pnts);
crdnts(1,:) = crvtr*angls_taus(1,:).*cos(angls_taus(1,:));
crdnts(2,:) = crvtr*angls_taus(1,:).*sin(angls_taus(1,:));
crdnts(3,:) = angls_taus(2,:);

% plot(angls_taus(1,:), angls_taus(2,:), '.', 'MarkerSize', 1); axis equal;
% 
% plot3(crdnts(1,:), crdnts(2,:), crdnts(3,:), '.', 'MarkerSize', .5)
% xlabel('x'); ylabel('y'); zlabel('z')
% axis square
% end