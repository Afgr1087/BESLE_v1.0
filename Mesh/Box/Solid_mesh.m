function Solid_mesh(a,b,c,nreg_a,nreg_b,nreg_c)

global POINTS_unit_cube POINTS ELEM_unit_cube ELEM FACES_unit_cube FACES
global NORMAL_VECTORS_unit_cube NORMAL_VECTORS nreg El_reg Volumes

space_x = a/nreg_a; % Space in x
space_y = b/nreg_b; % space in y
space_z = c/nreg_c; % space in z

% Scaling the unit cube
POINTS_unit_cube(:,2) = POINTS_unit_cube(:,2)*space_x; % x
POINTS_unit_cube(:,3) = POINTS_unit_cube(:,3)*space_y; % y
POINTS_unit_cube(:,4) = POINTS_unit_cube(:,4)*space_z; % z

points_aux = POINTS_unit_cube;
npoints = length(points_aux(:,1));
ind1 = 1;
ind2 = npoints;

elem_aux = ELEM_unit_cube;
nelem = length(elem_aux(:,1));
ind3 = 1;
ind4 = nelem;

face_aux = FACES_unit_cube;
nface = length(face_aux(:,1));
ind5 = 1;
ind6 = nface;

nv_aux = NORMAL_VECTORS_unit_cube;
nnv = length(nv_aux(:,1));
ind7 = 1;
ind8 = nnv;

cont = 0;
%% Volummes

Volumes = zeros(nreg,2);
Volumes(:,1) = [1:1:nreg];
Volumes(:,2) = space_x*space_y*space_z;

%% Creating the multiregions by levels
% Matrix of points and elements
for k=1:nreg_c % Third level in z axis 
    
    if (cont>0)
    	points_aux(:,4) = points_aux(:,4) + space_z;
    end

    for j=1:nreg_b % Second level in y axis
       
        if (j>1)
            points_aux(:,3) = points_aux(:,3) + space_y;
        end

        for i=1:nreg_a % First level in x axis
            
            if (cont==0) % only for the fisrt region
                % Matrix POINTS
                POINTS(ind1:ind2,:) = points_aux;
                ind1 = ind2 + 1;
                ind2 = ind2 + npoints;
                
                % Matrix ELEM
                ELEM(ind3:ind4,:) = elem_aux;
                ind3 = ind4 + 1;
                ind4 = ind4 + nelem;
            end                                        
            if (i==1 & cont>0) 
                
                % Matrix POINTS
                POINTS(ind1:ind2,:) = points_aux;
                ind1 = ind2 + 1;
                ind2 = ind2 + npoints;
                                
                % Matrix ELEM
                elem_aux(:,3:5) = elem_aux(:,3:5) + npoints;
                ELEM(ind3:ind4,:) = elem_aux;  
                ind3 = ind4 + 1;
                ind4 = ind4 + nelem;
            end
            if(i>1)
                
                % Matrix POINTS
                points_aux(:,2) = points_aux(:,2) + space_x;
                POINTS(ind1:ind2,:) = points_aux;
                ind1 = ind2 + 1;
                ind2 = ind2 + npoints;
                                
                % Matrix ELEM
                elem_aux(:,3:5) = elem_aux(:,3:5) + npoints;
                ELEM(ind3:ind4,:) = elem_aux;
                ind3 = ind4 + 1;
                ind4 = ind4 + nelem;
            end
            
            cont = cont + 1;
        
            % Matrix FACES
            FACES(ind5:ind6,:) = face_aux;
            ind5 = ind6 + 1;
            ind6 = ind6 + nface;
                
            % Matrix NORMAL_VECTORS
            NORMAL_VECTORS(ind7:ind8,:) = nv_aux;
            ind7 = ind8 + 1;
            ind8 = ind8 + nnv;
            
            % Matrix of Subregions
            El_reg(cont,:) = [cont,nelem];
                        
        end
        points_aux(:,2) = POINTS_unit_cube(:,2);
    end
    points_aux(:,3) = POINTS_unit_cube(:,3);
end

if (nreg>1)
    POINTS(:,1) = [1:npoints*nreg];
    ELEM(:,1) = [1:nelem*nreg];
    FACES(:,1) = [1:nface*nreg];
    NORMAL_VECTORS(:,1) = [1:nnv*nreg];
end
% FACE by element
ind=0;
for i=1:nface*nreg
    n_el = FACES(i,2);
    for j=1:n_el
        ind = ind + 1;
        ELEM(ind,2) = i;
    end
end

new_nv = [];
nv_aux = [];
ind2 = 0;
for i=1:length(FACES(:,1))
    nv_aux = NORMAL_VECTORS(i,:);
    nel_face = FACES(i,2);
    ind1 = ind2 + 1;
    ind2 = ind2 + nel_face;
    for j=ind1:ind2
        new_nv(j,:) = nv_aux;
    end
end

NORMAL_VECTORS = [];
NORMAL_VECTORS = new_nv;

vals = ELEM(:,1);
FACES = [];
FACES(:,1) = vals;
FACES(:,2) = 1;
ELEM(:,2) = vals;

end