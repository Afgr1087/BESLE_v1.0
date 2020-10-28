fprintf('3D MESH GENERATOR OF BOXES\n')
fprintf('Using 3-node triangular Discontinuous elements\n')
fprintf('By Andres F. Galvis\n')
fprintf('School of Mechanical Engineering\n')
fprintf('University of Campinas\n')
fprintf('Date: 19/02/2020\n')

global nreg

fprintf('\n 1. Triangle mesh generator\n')

% Total number of regions
nreg = nreg_a*nreg_b*nreg_c;
fprintf('\n 1.1 Number of regions: %d',nreg)

fprintf('\n\n 1.2 Generating the unit cube discretized \n')
% Unit cube discretized
Unit_cubic_mesh(nel_a,nel_b,nel_c);

fprintf('\n 1.3 Multiregion discretization \n')
Solid_mesh(a,b,c,nreg_a,nreg_b,nreg_c);

fprintf('\n2. Exporting data to fortran code\n')
Export_data;

