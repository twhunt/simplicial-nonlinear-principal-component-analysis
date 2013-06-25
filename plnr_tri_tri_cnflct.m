function cnflct = plnr_tri_tri_cnflct(...
    tri1_vrtx_crds_x, tri1_vrtx_crds_y, ...
    tri2_vrtx_crds_x, tri2_vrtx_crds_y, ...
    tri1_shrd_vertcs, tri2_shrd_vertcs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG
%hold off
%%subplot(2,1,2);
%plot(...
%   [tri1_vrtx_crds_x tri1_vrtx_crds_x(1)], ...
%   [tri1_vrtx_crds_y tri1_vrtx_crds_y(1)], '.-' , ...
%   [tri2_vrtx_crds_x tri2_vrtx_crds_x(1)], ...
%   [tri2_vrtx_crds_y tri2_vrtx_crds_y(1)], '*-')
% DEBUG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%if:
%tri1_shrd_vertcs(1) = 2 and tri2_shrd_vertcs(1) = 3
%then:
%vertex 2 of triangle 1 coincides with vertex 3 of triangle 2

numel_tri_shrd_vertcs = numel(tri1_shrd_vertcs);


if numel_tri_shrd_vertcs ~= numel(tri2_shrd_vertcs);
    error(['Arrays of shared triangle vertex indices' ...
        'must have the same number of elements']);
elseif numel_tri_shrd_vertcs > 3
    error(['Arrays of shared triangle vertices must have ' ...
        '0, 1, 2, or 3 elements']);
end

switch numel_tri_shrd_vertcs
    case 0
        
        
        tri1_vert_crds1 = ...
            [tri1_vrtx_crds_x(1); ...
             tri1_vrtx_crds_y(1)];
        
        tri1_vert_crds2 = ...
            [tri1_vrtx_crds_x(2); ...
             tri1_vrtx_crds_y(2)];

        tri1_vert_crds3 = ...
            [tri1_vrtx_crds_x(3); ...
             tri1_vrtx_crds_y(3)];

        
        tri2_vert_crds1 = ...
            [tri2_vrtx_crds_x(1); ...
             tri2_vrtx_crds_y(1)];
        
        tri2_vert_crds2 = ...
            [tri2_vrtx_crds_x(2); ...
             tri2_vrtx_crds_y(2)];

        tri2_vert_crds3 = ...
            [tri2_vrtx_crds_x(3); ...
             tri2_vrtx_crds_y(3)];

        
        cnflct = tri_tri_0_shrd_vrt_cnflct(...
            tri1_vert_crds1, tri1_vert_crds2, tri1_vert_crds3, ...
            tri2_vert_crds1, tri2_vert_crds2, tri2_vert_crds3);

    case 1
    
        %tri1_shrd_vertcs is a scalar
        
        
        %coordinates of the single vertex common to both triangles
        edge_tri_vert_crds = ...
            [tri1_vrtx_crds_x(tri1_shrd_vertcs); ...
             tri1_vrtx_crds_y(tri1_shrd_vertcs)]; 
        
        %treat triangle 1 as the candidate triangle in the sense of the
        %input arguments to tri_tri_1_shrd_vrt_cnflct
        switch tri1_shrd_vertcs
            case 1
                tri1_nonshrd_vrtcs = [2 3];
            case 2
                tri1_nonshrd_vrtcs = [1 3];
            case 3
                tri1_nonshrd_vrtcs = [1 2];
            otherwise
                error('Illegal entry in tri1_nonshrd_vrtcs')
        end
                
        edge_vert_crds = ...
            [tri1_vrtx_crds_x(tri1_nonshrd_vrtcs(1)); ...
             tri1_vrtx_crds_y(tri1_nonshrd_vrtcs(1))]; 

        cand_vert_crds = ...
            [tri1_vrtx_crds_x(tri1_nonshrd_vrtcs(2)); ...
             tri1_vrtx_crds_y(tri1_nonshrd_vrtcs(2))]; 

        %treat triangle 2 as the existing triangle in the sense of the
        %input arguments to tri_tri_1_shrd_vrt_cnflct
        switch tri2_shrd_vertcs
            case 1
                tri2_nonshrd_vrtcs = [2 3];
            case 2
                tri2_nonshrd_vrtcs = [1 3];
            case 3
                tri2_nonshrd_vrtcs = [1 2];
            otherwise
                error('Illegal entry in tri1_nonshrd_vrtcs')
        end
                
        tri_vert_crds2 = ...
            [tri2_vrtx_crds_x(tri2_nonshrd_vrtcs(1)); ...
             tri2_vrtx_crds_y(tri2_nonshrd_vrtcs(1))]; 

        tri_vert_crds3 = ...
            [tri2_vrtx_crds_x(tri2_nonshrd_vrtcs(2)); ...
             tri2_vrtx_crds_y(tri2_nonshrd_vrtcs(2))]; 
                
        cnflct = tri_tri_1_shrd_vrt_cnflct(...
            cand_vert_crds, edge_tri_vert_crds, edge_vert_crds, ...
            tri_vert_crds2, tri_vert_crds3);
        
    case 2
    
        %calling convention for tri_tri_2_shrd_vrt_cnflct is different than
        %tri_tri_0_shrd_vrt_cnflct and tri_tri_1_shrd_vrt_cnflct (for no
        %good reason
        %the nonshared vertex of one of the triangles is assumed to be the
        %origin, so translate all coordinates 
        
        %following if-elseif is equivalent to 
        %setdiff([1 2 3], tri1_shrd_vertcs)
        if tri1_shrd_vertcs(1) < tri1_shrd_vertcs(2)

            if tri1_shrd_vertcs(1) == 1 && tri1_shrd_vertcs(2) == 2
                tri1_nonshrd_vrtcs = 3;
            elseif tri1_shrd_vertcs(1) == 1 && tri1_shrd_vertcs(2) == 3
                tri1_nonshrd_vrtcs = 2;
            elseif tri1_shrd_vertcs(1) == 2 && tri1_shrd_vertcs(2) == 3
                tri1_nonshrd_vrtcs = 1;
            else
                error('Illegal entry in tri1_shrd_vertcs')
            end
            
        elseif tri1_shrd_vertcs(1) > tri1_shrd_vertcs(2)
            
            if tri1_shrd_vertcs(1) == 2 && tri1_shrd_vertcs(2) == 1
                tri1_nonshrd_vrtcs = 3;
            elseif tri1_shrd_vertcs(1) == 3 && tri1_shrd_vertcs(2) == 1
                tri1_nonshrd_vrtcs = 2;
            elseif tri1_shrd_vertcs(1) == 3 && tri1_shrd_vertcs(2) == 2
                tri1_nonshrd_vrtcs = 1;
            else
                error('Illegal entry in tri2_shrd_vertcs')
            end

        end
        
        %translate vertex coordinates
        %nonshared vertex of triangle gets translated to the origin
        tri1_vrtx_crds_x(tri1_shrd_vertcs) = ...
            tri1_vrtx_crds_x(tri1_shrd_vertcs) ...
            - tri1_vrtx_crds_x(tri1_nonshrd_vrtcs);
        
        tri1_vrtx_crds_y(tri1_shrd_vertcs) = ...
            tri1_vrtx_crds_y(tri1_shrd_vertcs) ...
            - tri1_vrtx_crds_y(tri1_nonshrd_vrtcs);
        
        %only need to translate the nonshared vertex of triangle 2
        if tri2_shrd_vertcs(1) < tri2_shrd_vertcs(2)

            if tri2_shrd_vertcs(1) == 1 && tri2_shrd_vertcs(2) == 2
                tri2_nonshrd_vrtcs = 3;
            elseif tri2_shrd_vertcs(1) == 1 && tri2_shrd_vertcs(2) == 3
                tri2_nonshrd_vrtcs = 2;
            elseif tri2_shrd_vertcs(1) == 2 && tri2_shrd_vertcs(2) == 3
                tri2_nonshrd_vrtcs = 1;
            else
                error('Illegal entry in tri2_shrd_vertcs')
            end
            
        elseif tri2_shrd_vertcs(1) > tri2_shrd_vertcs(2)
            
            if tri2_shrd_vertcs(1) == 2 && tri2_shrd_vertcs(2) == 1
                tri2_nonshrd_vrtcs = 3;
            elseif tri2_shrd_vertcs(1) == 3 && tri2_shrd_vertcs(2) == 1
                tri2_nonshrd_vrtcs = 2;
            elseif tri2_shrd_vertcs(1) == 3 && tri2_shrd_vertcs(2) == 2
                tri2_nonshrd_vrtcs = 1;
            else
                error('Illegal entry in tri2_shrd_vertcs')
            end

        end

        tri2_vrtx_crds_x(tri2_nonshrd_vrtcs) = ...
            tri2_vrtx_crds_x(tri2_nonshrd_vrtcs) ...
            - tri1_vrtx_crds_x(tri1_nonshrd_vrtcs);
        
        tri2_vrtx_crds_y(tri2_nonshrd_vrtcs) = ...
            tri2_vrtx_crds_y(tri2_nonshrd_vrtcs) ...
            - tri1_vrtx_crds_y(tri1_nonshrd_vrtcs);


        cnflct = plnr_tri_tri_2_shrd_vrt_cnflct(...
            tri1_vrtx_crds_x(tri1_shrd_vertcs), ...
            tri1_vrtx_crds_y(tri1_shrd_vertcs), ...
            tri2_vrtx_crds_x(tri2_nonshrd_vrtcs), ...
            tri2_vrtx_crds_y(tri2_nonshrd_vrtcs));
        
    case 3
    
        %triangles are identical
        cnflct = true;
end