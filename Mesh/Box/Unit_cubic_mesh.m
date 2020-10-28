%% Triangle_Mesh generator for unit cubic-rectangle geometries
function Unit_cubic_mesh(nel_a,nel_b,nel_c);
global FACES_unit_cube POINTS_unit_cube ELEM_unit_cube
global NORMAL_VECTORS_unit_cube BC_Face 

a = 1; % Dimension in x
b = 1; % Dimension in y
c = 1; % Dimension in z

% FACES = [face, N. elements]
FACES_unit_cube = [1 2*(nel_a*nel_c)
                   2 2*(nel_b*nel_c)
                   3 2*(nel_a*nel_c)
                   4 2*(nel_b*nel_c)
                   5 2*(nel_a*nel_b)
                   6 2*(nel_a*nel_b)];
     
% Each face is divided according to the number of elements

%% FACE 1
z1 = c;
z2 = c;
j=1;
k=1;
i=0;
space_z = c/nel_c; % Space between new points
space = a/nel_a; % space between new points
% Square divisions in the face
for jj=1:nel_c
    z2 = z2 - space_z; % x coordinate of the point 1
    POINTS_unit_cube(j:j+1,:) = [j,0,0,z1;j+1,0,0,z2];
    j=j+1;
    x1 = 0; % x coordinate of the point 2
    x2 = 0; % x coordinate of the point 2
    for ii=1:nel_a % Loop for new points
        j=j+1;
        x1 = x1 + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x1,0,z1];
        j = j+1;
        x2 = x2 + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x2,0,z2]; 
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k = k+2;
    j=j+1;
    z1 = z1 - space_z; 
end
N_points_face(1) = j-1; % Number of points in FACE 1
% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% FACE 2

z1 = c;
z2 = c;

space_z = c/nel_c; % Space between new points
space = b/nel_b; % space between new points
% Square divisions in the face
for jj=1:nel_c
    z2 = z2 - space_z; % x coordinate of the point 1
    POINTS_unit_cube(j:j+1,:) = [j,a,0,z1;j+1,a,0,z2];
    j=j+1;
    ya = 0; % x coordinate of the point 2
    yb = 0; % x coordinate of the point 2
    for ii=1:nel_b % Loop for new points
        j=j+1;
        ya = ya + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,a,ya,z1];
        j = j+1;
        yb = yb + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,a,yb,z2]; 
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k=k+2;
    j=j+1;
    z1 = z1 - space_z; 
end
N_points_face(2) = j-1; % Number of points in FACE 2

% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% FACE 3
% We need to copy all POINTS in the FACE 1 to the FACE 3, but
% the y coordinate plus b

pos_1 = N_points_face(1); % Number of points in FACE 3 the same of FACE 1
points_aux = POINTS_unit_cube(1:pos_1,:);

pos_2 = N_points_face(2);
points_aux(:,1) = [pos_2+1:pos_2+pos_1]; % New numeration
points_aux(:,3) = points_aux(:,3) + b; % New coordinate in y axis

POINTS_unit_cube(pos_2+1:pos_2+pos_1,:) = points_aux; % Points in FACE 3
N_points_face(3) = pos_2+pos_1; % Number of points in FACE 2
% quads
for jj=1:nel_c
    for ii=1:nel_a % Loop for new points
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k=k+2;
end

% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% FACE 4
% We need to copy all POINTS in the FACE 2 to the FACE 4, but
% the x coordinate minus a

pos_1 = N_points_face(1); % Number of points in FACE 1
pos_2 = N_points_face(2); % Number of points in FACE 4 the same of FACE 2
points_aux = POINTS_unit_cube(pos_1+1:pos_2,:);

pos_3 = N_points_face(3); % Number of points in FACE 3
points_aux(:,1) = [pos_3+1:pos_3+(pos_2-pos_1)]; % New numeration
points_aux(:,2) = points_aux(:,2) - a; % New coordinate in x axis

POINTS_unit_cube(pos_3+1:pos_3+(pos_2-pos_1),:) = points_aux; % Points in FACE 4

% quads
for jj=1:nel_c
    for ii=1:nel_b % Loop for new points
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k=k+2;
end

% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% FACE 5

j = length(POINTS_unit_cube(:,1)) + 1;

x1 = 0;
x2 = 0;

space_x = a/nel_a; % Space between new points
space = b/nel_b; % space between new points
for jj=1:nel_a
    x2 = x2 + space_x; % x coordinate of the point 1
    POINTS_unit_cube(j:j+1,:) = [j,x1,0,0;j+1,x2,0,0];
    j=j+1;
    ya = 0; % x coordinate of the point 2
    yb = 0; % x coordinate of the point 2
    for ii=1:nel_b % Loop for new points
        j=j+1;
        ya = ya + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x1,ya,0];
        j = j+1;
        yb = yb + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x2,yb,0]; 
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k=k+2;
    j=j+1;
    x1 = x1 + space_x; 
end
% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% FACE 6

x1 = 0;
x2 = 0;

space_x = a/nel_a; % Space between new points
space = b/nel_b; % space between new points
for jj=1:nel_a
    x2 = x2 + space_x; % x coordinate of the point 1
    POINTS_unit_cube(j:j+1,:) = [j,x1,0,c;j+1,x2,0,c];
    j=j+1;
    ya = 0; % x coordinate of the point 2
    yb = 0; % x coordinate of the point 2
    for ii=1:nel_b % Loop for new points
        j=j+1;
        ya = ya + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x1,ya,c];
        j = j+1;
        yb = yb + space; % New x coordinate
        POINTS_unit_cube(j,:) = [j,x2,yb,c]; 
        % FACES
        i=i+1;
        QUADS(i,:) = [i  k  k+1  k+3  k+2];
        k = k+2;
    end
    k=k+2;
    j=j+1;
    x1 = x1 + space_x; 
end
% Triangle divisions in each square
ntri = size(QUADS,1);
ind=0;
node = [2,4];
for ii=1:ntri
   node1 = QUADS(ii,3); node2 = QUADS(ii,5); 
   for jj = 1:2
       ind = ind + 1;
       kk = node(jj);
       node3 = QUADS(ii,kk);
       ELEM_unit_cube(ind,1) = ind;
       ELEM_unit_cube(ind,3:5) = [node1, node2, node3];
   end
end

%% ----------------------------------------------------------------------- %
% NORMAL_VECTORS : [Face, nx, ny, nz]                            %
% ----------------------------------------------------------------------- %

NORMAL_VECTORS_unit_cube = [1  0 -1  0
                            2  1  0  0
                            3  0  1  0
                            4 -1  0  0
                            5  0  0 -1
                            6  0  0  1];

end