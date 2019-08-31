function [U, res, sol] = BoundUSys(bound,f,UODESys_input,magnitude,d,epsilon,sigma,verbose,SOSPsolver,SDPsolver)
% Compute an upper or lower bound on the specifid polynoamial magnitude in a N-dimesional hydrodinamic-type system
% after reduction to a n-dimensional uncertain system with (n<N)

% bound         'U' to compute upper bound, 'L' to compute lower bound
% f             File containing N_ijk, L_ij and B_i            
% UODESys_input {n,[g1,g2,g3],output_file} inputs for UODEsys
% magnitude     Magnitude to be bouunded
% d             Half the degree of the auxiliary function
% epsilon       Noise intensity
% sigma
% verbose       0 or 1
% SOSPsolver    'yalmip' and 'spotless' are available
% SDPsolver     SDP solver. Only sedumi and mosek are implemented for spotless

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics


% Load N_ijk, L_ij and B_i 
load(f,'N_ijk','L_ij','B_i');
N = length(B_i);
n = UODESys_input{1};
g = UODESys_input{2};
UODESys_file = UODESys_input{3};
% Assign default arguments
if nargin < 6, epsilon = 0; end
if nargin < 7, sigma = eye(n); end
if nargin < 8, verbose = 0; end
if nargin < 9, SOSPsolver = 'yalmip'; end
if nargin < 10,SDPsolver = 'mosek'; end

% Reduced the initial finite dimensional system to an uncertain system if not done yet
if not(isfile(UODESys_file)) UODESys(n,g,f,UODESys_file,SDPsolver,verbose); end
load(UODESys_file);
% Check that c_1, c_2 and c_3 obtained with UODESys are >= 0
if (c_1_solution < 0 | c_2_solution < 0 | c_3_solution < 0)
  disp("ERROR: The input file to BoundUncertain must satisfy c_1_solution >= 0, c_2_solution >= 0 and c_3_solution >= 0.");
  U = inf;
  res = inf;
  return
end
n = length(B_hat_x_final);
disp("BoundUSys: Finding a bound for the "+n+"-dimensional uncertain system...")

switch SOSPsolver
    case 'yalmip'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%% YALMIP %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("BoundUSys: Using YALMIP.");

clear('yalmip')
% Define system independet variables
a = sdpvar(n,1);
q  = sdpvar(1);
q2 = q^2;
z = sdpvar(n,1);
z0 = sdpvar(1);
% Define system dynamics
if (N_hat_xyz_final == 0)
    f = L_hat_xy_final*a + B_hat_x_final;
else
    f = Naa(N_hat_xyz_final,a) + L_hat_xy_final*a + B_hat_x_final;
end

% Define average of interest applying proper transformation
phi = magnitude(eigenvectors(1:n,1:n)*a,q2);

% Define auxiliary functional
sdpvar c;
[p,pc] = polynomial([a;q2],2*d-1);
V = c*(a'*a/2 + q2)^d + p;

% Define upper bound for ||\Theta||^2
P = c_1_solution*q2 + 2*c_2_solution*(a'*a)*q2 + 4*c_3_solution*q2^2;

% Define D
sdpvar U;
if epsilon == 0
    % No noise. Deterministic case
    D = -(P*z0^2 + z'*z)*(jacobian(V,a)*f + jacobian(V,q)*.5*KAPPA*q + phi - U) - 2*P*z0*(jacobian(V,a) - jacobian(V,q)*(0.5/q)*a')*z;
elseif size(sigma,1) == n   
    % Noise is added directly to the uncertain system
    % Compute matrix S for the system with state vector (b_1,b_2,...,b_n)
    S = sigma*sigma';
    D = -(P*z0^2 + z'*z)*(jacobian(V,a)*f + jacobian(V,q)*.5*KAPPA*q + phi + epsilon*div( S*jacobian(V,a)',a) - U) - 2*P*z0*(jacobian(V,a) - jacobian(V,q)*(0.5/q)*a')*z;
else
    disp("BoundUsys: ERROR. Not valid dimension for sigma. Sigma must be a nxn matrix, where n is the dimension of the uncertain system.")
    return
end

% SDP solver options
ops = sdpsettings('solver',SDPsolver,'verbose',verbose);
ops.mosek.MSK_IPAR_NUM_THREADS = 2;
ops.savesolveroutput = 1;
ops.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;
ops.sos.scale = 1;

switch bound
    case 'U'
        % Find an Upper bound
        constraints = [sos(jacobian(V,q)/q); sos(D)];
        [sol,~,~,res] = solvesos(constraints,U,ops,[U;pc;c]);
    case 'L'
        % Find a Lower bound
        constraints = [sos(jacobian(V,q)/q); sos(-D)];
        [sol,~,~,res] = solvesos(constraints,-U,ops,[U;pc;c]);
    otherwise
        disp("BoundUSys: ERROR. First argument not valid.");
        disp("bound = 'U' for upper bound and bound = 'L' for lower bound.");
        return
end 

U = value(U);
res = norm(res,Inf);

case 'spotless'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%% SPOTLess %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("BoundUSys: Using SPOTLess.");

% SPOTLess options
if strcmp(SDPsolver,'mosek')
    solver = @spot_mosek;
elseif strcmp(SDPsolver,'sedumi')
    solver = @spot_sedumi;
else
   disp("BoundUSys: ERROR. Not valid SDP solver."); 
end    
opts = spot_sdp_default_options(); % Create structure
opts.verbose = verbose; % Verbosity for SDP solver
opts.solver_options.mosek.MSK_IPAR_NUM_THREADS = 2; % Set number of processors to 2
opts.solver_options.mosek.MSK_IPAR_INTPNT_MULTI_THREAD = 'MSK_ON';  % enable parallelization in mosek
opts.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;

% Initialise SoS program
prog = spotsosprog;
a  = msspoly('a',n);               % Define state vector
q  = msspoly('q',1);
q2 = q^2;
z0  = msspoly('zo',1);
z  = msspoly('z',n);
prog = prog.withIndeterminate([a;q;z0;z]); % Define independent variables

% Define system dynamics
linear_term = L_hat_xy_final*a;
if (N_hat_xyz_final == 0)
    f = linear_term;
else
    f = NaaSPOT(N_hat_xyz_final,a) + linear_term + B_hat_x_final;
end

% Define average of interest applying proper transformation
phi = magnitude(eigenvectors(1:n,1:n)*a,q2);

% Create V: exponents of monomials, monomials, coefficients and V
Vpow = my_monpowers(n+1,2*d-1);
Vpow = [Vpow(:,1:end-1), 2*Vpow(:,end)];
np   = size(Vpow,1);
Vm = recomp([a;q], Vpow, speye(size(Vpow,1)) );
[prog,Vc] = prog.newFree(size(Vpow,1));
[prog,c] = prog.newFree(1); % coefficient of energy^d
p = Vc'*Vm;
energy = a'*a/2 + q2;
E = c*energy^d;
V = E + p;

Vpow(:,end) = Vpow(:,end)/2; % get powers of q^2 instead of powers of q
monomials_to_diff = Vpow(:,end)>0;
coefficients_of_dV_by_dq2 = Vpow(monomials_to_diff,end).*Vc(monomials_to_diff);
powers_in_dV_by_dq2_variables_a_and_q  = [Vpow(monomials_to_diff,1:end-1), 2*(Vpow(monomials_to_diff,end)-1)];
monomials_in_dV_by_dq2 = recomp([a; q],powers_in_dV_by_dq2_variables_a_and_q, speye(size(powers_in_dV_by_dq2_variables_a_and_q,1)));
dV_by_dq2 = c*d*(a'*a/2+q2)^(d-1) + coefficients_of_dV_by_dq2'*monomials_in_dV_by_dq2;

% Variable to be bounded
[prog,U] = prog.newFree(1);

% SOS constraint
P = c_1_solution*q2 + 2*c_2_solution*(a'*a)*q2 + 4*c_3_solution*q2^2;
if epsilon == 0
    D = -(P*z0^2 + z'*z)*(diff(V,a)*f + dV_by_dq2*q2*KAPPA + phi - U) - 2*P*z0*(diff(V,a) - dV_by_dq2*a')*z  
elseif size(sigma,1) == n   
    % Noise is added directly to the uncertain system
    % Compute matrix S for the system with state vector (b_1,b_2,...,b_n)
    S = sigma*sigma';
    D = -(P*z0^2 + z'*z)*(diff(V,a)*f + dV_by_dq2*KAPPA*q2 + phi + epsilon*divSPOT( S*diff(V,a)',a) - U) - 2*P*z0*(diff(V,a) - dV_by_dq2*a')*z;
else
    disp("BoundUsys: ERROR. Not valid dimension for sigma. Sigma must be a nxn matrix, where n is the dimension of the uncertain system.")
    return
end

% Candidate monomials for SOS decomposition
% Exploit symmetry: 2 blocks
powers = my_monpowers(n+1,d);
% powers = [powers(:,1:end-1), 2*powers(:,end)];
h = recomp([a;q], powers, speye(size(powers,1)) );
monoms = {h};  

% Minimise upper bound
switch bound
    case "U"
        % Set constraints (Upper bound)
        prog = prog.withSOS(D);
        prog = prog.withSOS(dV_by_dq2);
   %     prog = prog.withSparseSOS(dV_by_dq2,monoms);     
        sol = prog.minimize(U,solver,opts);
        % Compute D(a) - h^TQh
        res1 =   D - sol.gramMonomials{1}'*double(sol.eval(sol.gramMatrices{1}))*sol.gramMonomials{1};
    case "L"
        % Set constraints (Lower bound)
        prog = prog.withSOS(-D);
        prog = prog.withSOS(dV_by_dq2);
        sol = prog.minimize(-U,solver,opts); 
        % Compute D(a) - h^TQh
        res1 =   -D - sol.gramMonomials{1}'*double(sol.eval(sol.gramMatrices{1}))*sol.gramMonomials{1};      
    otherwise
        disp("BoundUSys ERROR: First argument not valid.");
        disp("bound = 'U' for upper bound and bound = 'L' for lower bound.");
        return
        
end

U   = double(sol.eval(U)); 
% Compute residual
res1 = sol.eval(res1);
res1 = norm(res1.coeff,Inf);
res2 =   dV_by_dq2 - sol.gramMonomials{2}'*double(sol.eval(sol.gramMatrices{2}))*sol.gramMonomials{2};
res2 = sol.eval(res2);
res2 = norm(res2.coeff,Inf);
res = norm([res1,res2],Inf);
    end

