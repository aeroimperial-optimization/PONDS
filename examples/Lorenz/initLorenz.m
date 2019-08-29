function initLorenz(beta,sigma,rho,scaling_factor,file)
% Compute N_ijk, L_ij and B_i for Lorenz atractor
%
% beta
% sigma
% rho
% scaling_factor    Scaling factor
% file              Destination file

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

N = 3; % Lorenz problem is 3D

% Build N_ijk
N_ijk = zeros(N,N,N);
N_ijk(2,3,1) = -.5;
N_ijk(2,1,3) = -.5;
N_ijk(3,2,1) =  .5;
N_ijk(3,1,2) =  .5;

% Build L_ijk
L_ij = [-sigma sigma 0; rho -1 0; 0 0 -beta];

% There are no constant terms
B_i = zeros(N,1);

% Apply rescaling 
N_ijk = N_ijk*scaling_factor;

save(file,"N_ijk","L_ij","B_i");
