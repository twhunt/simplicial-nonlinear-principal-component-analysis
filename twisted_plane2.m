function [crdntsND] = twisted_plane2(crdnts2D, max_twist_angls)
%Twists a 2D plane into N dimensions, where N=numel(max_twist_angls) + 2
%To twist a 2D plane into a 3D plane, use the map
%[x; y] ->[x; y*cos(angl(x)); y*sin(angl(x))],
%where angl(x) = (x/x_max)*angl_max, that is treat [x y] as a
%coordinate vector for the basis {[1;0;0], [0; cos; sin]}
%In higher dimensions, the second basis vector is a point on the
%unit hypersphere

for k=2:-1:1
    max_ordnt(k,1) = max(crdnts2D(k, :));
end

num_pnts = size(crdnts2D, 2);
num_twist_angls = numel(max_twist_angls);

out_dmnsn = 2 + num_twist_angls;
crdntsND = zeros(out_dmnsn, num_pnts);

twist_fctr = max_twist_angls/max_ordnt(1);


%Basis at current point is I(:,1) and second column of
crdntsND(1, :) = crdnts2D(1, :);
for pnt_i=1:num_pnts
    
    tmp_angls = twist_fctr*crdnts2D(1, pnt_i);
    unit_vctr = unit_vector(tmp_angls);        
    crdntsND(2:end, pnt_i) = crdnts2D(2, pnt_i)*unit_vctr;
    
end

%[dbg_U dbg_S dbg_V] = svd(dbg_hshldr_mtrx_clmn);

% plot3(...
%     crdntsND(1, :), ...
%     crdntsND(2, :), ...    
%     crdntsND(3, :), ...
%     '.')
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = unit_vector(angls)
%returns a vector on the surface of the unit hypersphere
%hypersphere is parameterized by angles as follows
%R^2: u2 = [cos(angl(1)); sin(angl(1))]
%R^3: u3 = [cos(angl(2))*u2; sin(angl(2)]
%R^4: u4 = [cos(angl(3))*u3; sin(angl(3)]
%etc.
%Note the difference between standard spherical coordinates in R^3, where
%the third coordinate is the angle made with the z axis (and not the xy
%plane)

num_angls = numel(angls);
if num_angls == 0
    u = zeros(num_angls, 1);
else
    
    dmnsn = num_angls+1;
    u = zeros(dmnsn, 1);
    u(1) = cos(angls(1));
    u(2) = sin(angls(1));

    for crdnt_i=3:dmnsn
        
        u(crdnt_i) = sin(angls(crdnt_i-1));
        u(1:(crdnt_i-1)) = cos(angls(crdnt_i-1))*u(1:(crdnt_i-1));
        
    end
    
end

end