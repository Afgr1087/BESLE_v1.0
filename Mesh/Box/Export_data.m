%% ==========================[  EXPORT DATA  ]========================== %%

global FACES POINTS ELEM NORMAL_VECTORS BC_Face nreg El_reg Volumes

filename = strcat(file,'.dat');

[f p]=uiputfile(filename,'Save as Input_data.DAT');
fid = fopen(f, 'wt');
format long

fprintf(fid,'%10.0f\n',nreg);

% m = length(FACES(:,1));
% n = length(FACES(1,:));
% fprintf(fid,'%10.0f',m);
% fprintf(fid,'%10.0f\n',n);
% 
% for i = 1:size(FACES,1)
%     fprintf(fid,'%10.0f',FACES(i,:));
%     fprintf(fid,'\n');
% end

m = length(POINTS(:,1));
n = length(POINTS(1,2:4));
fprintf(fid,'%10.0f',m);
fprintf(fid,'%10.0f\n',n);

for i = 1:size(POINTS,1)
    fprintf(fid,'%30.15f',POINTS(i,2:4));
    fprintf(fid,'\n');
end

m = length(ELEM(:,1));
n = length(ELEM(1,3:5));
fprintf(fid,'%10.0f',m);
fprintf(fid,'%10.0f\n',n);

for i = 1:size(ELEM,1)
    fprintf(fid,'%10.0f',ELEM(i,3:5));
    fprintf(fid,'\n');
end

m = length(NORMAL_VECTORS(:,1));
n = length(NORMAL_VECTORS(1,2:4));
fprintf(fid,'%10.0f',m);
fprintf(fid,'%10.0f\n',n);

for i = 1:size(NORMAL_VECTORS,1)
    fprintf(fid,'%30.15f',NORMAL_VECTORS(i,2:4));
    fprintf(fid,'\n');
end

m = length(El_reg(:,1));
n = length(El_reg(1,2));
fprintf(fid,'%10.0f',m);
fprintf(fid,'%10.0f\n',n);

for i = 1:size(El_reg,1)
    fprintf(fid,'%10.0f',El_reg(i,2));
    fprintf(fid,'\n');
end

% m = length(Volumes(:,1));
% n = length(Volumes(1,:));
% fprintf(fid,'%10.0f',m);
% fprintf(fid,'%10.0f\n',n);
% 
% for i = 1:nreg
%     fprintf(fid,'%10.15f             ',Volumes(i,2));
%     fprintf(fid,'\n');
% end

fprintf(fid,'end\n');
fclose(fid);

%% Export *,vtk

filename2 = strcat(file,'.vtk');

[f p]=uiputfile(filename2,'Save as Input_data.DAT');
fid = fopen(f, 'wt');
format long

fprintf(fid,'%s\n','# vtk DataFile Version 2.3');
fprintf(fid,'%s\n',filename2);
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n',' ');
fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');

m = length(POINTS(:,1));
fprintf(fid,'POINTS%20d double\n',m);

for i = 1:size(POINTS,1)
    fprintf(fid,'%30.15f',POINTS(i,2:4));
    fprintf(fid,'\n');
end

m = length(ELEM(:,1));
fprintf(fid,'CELLS%20d%20d\n',m,m*4);

for i = 1:size(ELEM,1)
    fprintf(fid,'%10.0f',[3 ELEM(i,3:5)-1]);
    fprintf(fid,'\n');
end

fprintf(fid,'CELL_TYPES%20d\n',m);

for i = 1:size(ELEM,1)
    fprintf(fid,'%10d',5);
    fprintf(fid,'\n');
end