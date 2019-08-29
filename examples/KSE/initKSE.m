function [] = initKSE(L,N,scaling_factor,file)
% Compute N_ijk, L_ij and B_i for N-dimensional truncated 
% Kuramoto-Sivashinsky equation
%
% L                 Lenght scale in KSE
% N                 Number of modes in galkerkin truncation
% scaling_factor    Scaling factor
% file              Destination file

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

% Build N_ijk
N_ijk = zeros(N,N,N);
for i = 1:N
    for m = 1:N-i
       N_ijk(i,m,m+i) = N_ijk(i,m,m+i) + .5/sqrt(pi*L)*(i/L);
    end
    for m = 1:i-1
       N_ijk(i,m,i-m) = N_ijk(i,m,i-m) + -.25/sqrt(pi*L)*(i/L);
    end
end

% Build L_ijk
L_ij = zeros(N,N);
for i = 1:N
    L_ij(i,i) = (i/L)^2*(1-(i/L)^2);
end

% Build B_i
B_i = zeros(N,1);

% Apply rescaling 
N_ijk = N_ijk*scaling_factor;

save(file,"N_ijk","L_ij","B_i");
