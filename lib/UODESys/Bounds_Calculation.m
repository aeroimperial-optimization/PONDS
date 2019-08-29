function [] = Bounds_Calculation(n,g,output_file,solver,verbose)
% Script to calculate coefficient matrix for CHI and bounds for GAMMA and THETA.

% Written by Mayur Venkatram Lakshmi (May 2017)
% Imperial College London - Department of Aeronautics

% Important Note: 
% GAMMA and THETA in this MATLAB script represent GAMMA_hat and 
% THETA_hat, which are GAMMA and THETA after coordinate transformation.
% Similarly, CHI in this MATLAB script represents CHI_hat.

%% Reading in variables produced when running Variable_Transformation.m.
load('transformed_arrays.mat');
g
g1 = g(1); g2 = g(2); g3 = g(3);

%% Defining vectors for indexing.
u_index = 1:n;
v_index = n+1:N; 

%% Truncation of transformed arrays to obtain new uncertain dynamical system.

N_hat_xyz_final = N_hat_xyz(u_index,u_index,u_index);
L_hat_xy_final = L_hat_xy(u_index, u_index);
B_hat_x_final = B_hat_x(u_index);

disp('N_hat_xyz_final (N_ijk after variable transformation and truncation): ')
disp(N_hat_xyz_final)

disp('L_hat_xy_final (L_ij after variable transformation and truncation): ')
disp(L_hat_xy_final)

disp('B_hat_x_final (B_i after variable transformation and truncation): ')
disp(B_hat_x_final)

%% Determining coefficient matrix for CHI and a bound for GAMMA.
% L_hat_xy_sym is symmetric and diagonal from of L_hat_xy.
L_hat_xy_sym = 0.5*(L_hat_xy+ L_hat_xy');

% L_hat_xy_prime is the coefficient matrix for the GAMMA term (after
% coordinate transformation).
L_hat_xy_prime = L_hat_xy_sym(v_index,v_index);
disp('L_hat_xy_prime (The coefficient matrix of GAMMA_hat): ' )
disp(L_hat_xy_prime)

% The following two matrices are used to determine the coefficient matrix
% for the CHI term:
L_hat_xy_2_prime = L_hat_xy_sym(u_index,v_index);
L_hat_xy_3_prime = L_hat_xy_sym(v_index,u_index);

% L_hat_xy_star is the coefficient matrix for the CHI term (after
% coordinate transformation).
L_hat_xy_star = L_hat_xy_2_prime + L_hat_xy_3_prime';

% lambda_{n+1} - The (n+1)th eigenvalue of matrix A = 2*(L_ij + L_ji)
% It is used in determining the bound for GAMMA.
eigenvalue_n_plus_1 = eigenvalues(n+1,n+1);
disp('(n+1)th eigenvalue of A = 2*(L_ij + L_ji): ')
disp(eigenvalue_n_plus_1)

% Determining value of KAPPA, where GAMMA_hat <= KAPPA * (q_hat)^2.
KAPPA = 0.5*eigenvalues(n+1,n+1);
disp('KAPPA: ')
disp(KAPPA)

%% Determing bound for ||THETA|| using SOS methods.

% Recall that THETA represents THETA_hat.

% Coefficient matrix for quadratic term in Theta_x.
N_hat_xyz_tilde_array = cell(1,n);
for x = 1:n
    N_hat_xyz_tilde_array{1,x} = N_hat_xyz_array{1,x}(v_index,v_index);
end

% Coefficient matrix for bilinear term in Theta_x.
N_hat_xyz_2_tilde_array = cell(1,n);
for x = 1:n
    N_hat_xyz_2_tilde_array{1,x} = ...
        N_hat_xyz_array{1,x}(v_index,u_index)+ ...
        N_hat_xzy_array{1,x}(v_index,u_index);
end

% Coefficient matrix for linear term in Theta_x.
L_hat_xy_tilde_array = cell(1,n);
for x = 1:n
    L_hat_xy_tilde_array{1,x} = L_hat_xy_array{1,x}(v_index);
end

% Defining vectors of symbolic variables using YALMIP.
% NOTE: u and v are NEW (transformed) variables ie: 
% yalmip u = u_hat and yalmip v = v_hat.
u = sdpvar(1,n);
v = sdpvar(1,N-n);

% SUM (N_xyz_tilde * v_y * v_z) -- QUADRATIC TERM in Theta_x.
% Row cell array - Each cell corresponds to an x index value.
P1_array = cell(1,n);
for x = 1:n
    P1_array{x} = 0;
    for y = 1:N-n
        for z = 1:N-n
            P1_array{x} = P1_array{x} + ...
                N_hat_xyz_tilde_array{x}(y,z)*v(y)*v(z);
        end
    end
end

% SUM (N_xyz_2_tilde * u_z * v_y) -- BILINEAR TERM in Theta_x.
% Row cell array - Each cell corresponds to an x index value.
P2_array = cell(1,n);
for x = 1:n
    P2_array{x} = 0;
    for y = 1:N-n
        for z = 1:n 
            P2_array{x} = P2_array{x} + ...
                N_hat_xyz_2_tilde_array{x}(y,z)*u(z)*v(y);
        end
    end
end

% SUM (L_xy_tilde * v_y) -- LINEAR TERM in Theta_x.
% Row cell array - Each cell corresponds to an x index value.
P3_array = cell(1,n);
for x = 1:n
    P3_array{x} = 0;
    for y = 1:N-n
        P3_array{x} = P3_array{x} + L_hat_xy_tilde_array{x}(y)*v(y);
    end
end

% THETA = P_1 + P_2 + P_3
% NOTE: This YALMIP THETA is actually THETA_hat in derivation. ie: It is
% the THETA in new (transformed) variables.
Theta = cell(1,n);
for x = 1:n
    Theta{x} = P1_array{x} + P2_array{x} + P3_array{x};
end

% Another sdpvar variable -- (norm_Theta_sq) = ||THETA||^2.
norm_Theta_sq = 0;
for x = 1:n
    norm_Theta_sq = norm_Theta_sq + (Theta{x})^2;
end

% Another sdpvar variable -- norm_u_sq = || u ||^2.
norm_u_sq = 0;
for x = 1:length(u)
    norm_u_sq = norm_u_sq + (u(x))^2;
end

% Another sdpvar variable -- norm_v_sq = || v ||^2.
norm_v_sq = 0;
for x = 1:length(v)
    norm_v_sq = norm_v_sq + (v(x))^2;
end

% Defining symbolic decision variables using YALMIP. We want to find 
% optimum values of c_1, c_2 and c_3.
sdpvar c_1 c_2 c_3  

% Polynomial bound.
P = 0.5*c_1*norm_v_sq + c_2*norm_u_sq*norm_v_sq + c_3*(norm_v_sq)^2; 

% Setting up sum-of-squares (SOS) constraint: P - ||THETA||^2 = SOS.
Constraint = [sos(P - norm_Theta_sq),c_1>=1e-6,c_2>=1e-6,c_3>=1e-6];    

% We want to minimise the following objective function:
Objective = g1*c_1 + g2*c_2 + g3*c_3;

% Solve the optimization problem using YALMIP.

% Z (upper case) is a vector of monomials - NOT to be confused with index z
% P - ||THETA||^2 = Z^T * Q * Z.
ops = sdpsettings('solver',solver,'verbose',verbose);
ops.mosek.MSK_IPAR_NUM_THREADS = 6;
ops.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;
[sol,Z,Q,res] = solvesos(Constraint,Objective,ops,[c_1,c_2,c_3]); 

% Store the optimal values of c_1, c_2 and c_3 in new variables.
c_1_solution = value(c_1);
c_2_solution = value(c_2);
c_3_solution = value(c_3);

% Display the values of c_1, c_2 and c_3 in the command window.
disp('c1: ')
disp(c_1_solution)
disp('c2: ')
disp(c_2_solution)
disp('c3: ')
disp(c_3_solution)
disp('residual: ')
disp(res)

% Eigenvalues of Q. Has as many elements as the number of monomials in Z
Q_eigenvalues = eig(Q{1});

%% Save important variables in a file OUTPUT.mat.
save(output_file,'c_1_solution','c_2_solution','c_3_solution',...
    'L_hat_xy_star','KAPPA','L_hat_xy_prime','eigenvectors',...
    'eigenvalues','L_hat_xy','N_hat_xyz','B_hat_x','Q_eigenvalues', ...
    'N_hat_xyz_final','L_hat_xy_final','B_hat_x_final')
