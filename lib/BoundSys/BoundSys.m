function [U,res,sol] = BoundSys(bound,f,magnitude,d,epsilon,sigma,verbose,symmetries,SOSPsolver,SDPsolver)
% Compute an upper or lower bound on the specified polynomial magnitude in a N-dimesional hidrodynamic-type system
%
% bound         'U' to compute upper bound, 'L' to compute lower bound
% f             File containing N_ijk, L_ij and B_i
% magnitude     Magnitude to be bounded
% d             Half the degree of the auxiliary function
% epsilon       Noise intensity
% sigma
% verbose       0 or 1
% symmetries    Row vector containing the diag(\Lambda)
% SOSPsolver    'yalmip' and 'spotless' are available
% SDPsolver     SDP solver. Only sedumi and mosek are implemented for spotless

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

% Load N_ijk, L_ij and B_i
load(f,'N_ijk','L_ij','B_i');
N = length(B_i);

% Assign default arguments
if nargin < 5, epsilon = 0; end
if nargin < 6, sigma = eye(N); end
if nargin < 7, verbose = 0; end
if nargin < 8, symmetries = 0; end
if nargin < 9, SOSPsolver = 'yalmip'; end
if nargin < 10,SDPsolver = 'mosek'; end

disp("BoundSys: Finding a bound for the "+N+"-dimensional stochastic system...");

switch SOSPsolver
    case 'yalmip'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%% YALMIP %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("BoundSys: Using YALMIP");

clear('yalmip')
a = sdpvar(N,1); % Define state vector
% Define system dynamics
if (N_ijk == 0)
    f = L_ij*a + B_i;
else
    f = Naa(N_ijk,a) + L_ij*a + B_i;
end

% Define average of interest
phi = magnitude(a);

% Define auxiliary fuctional
if  abs(symmetries) == 1
        disp("BoundSys: Valid symmetry.")
        [p,pc] = symmpoly(a,2*d-1,symmetries);
else
        [p,pc] = polynomial(a,2*d-1);
end    
sdpvar c;
V = c*(a'*a)^d + p;

% SDP solver options
ops = sdpsettings('solver',SDPsolver,'verbose',verbose);
ops.mosek.MSK_IPAR_NUM_THREADS = 2;
ops.savesolveroutput = 1;
ops.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;

sdpvar U;
if epsilon == 0
    D = U - jacobian(V,a)*f - phi;
else
    D = U - phi - jacobian(V,a)*f - epsilon*div(sigma*sigma'*(jacobian(V,a))',a);
end

switch bound
    case "U"
        % Find an Upper bound
        [sol,~,~,res] = solvesos(sos(D),U,ops,[U;pc;c]);

    case "L"
        % Find a Lower bound
        [sol,~,~,res] = solvesos(sos(-D),-U,ops,[U;pc;c]);
            
    otherwise
        disp("BoundSys: ERROR. First argument not valid.");
        disp("bound = 'U' for upper bound and bound = 'L' for lower bound.");
        return
end

U = value(U);
res = norm(res);
case 'spotless'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%% SPOTLess %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("BoundSys: Using SPOTLess");

% SPOTLess options
if strcmp(SDPsolver,'mosek')
    solver = @spot_mosek;
elseif strcmp(SDPsolver,'sedumi')
    solver = @spot_sedumi;
else
   disp("BoundSys: ERROR. Not valid SDP solver."); 
end    
opts = spot_sdp_default_options(); % Create structure
opts.verbose = verbose; % Verbosity for SDP solver
opts.solver_options.mosek.MSK_IPAR_NUM_THREADS = 2; % Set number of processors to 2
opts.solver_options.mosek.MSK_IPAR_INTPNT_MULTI_THREAD = 'MSK_ON';  % enable parallelization in mosek
opts.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;

% Initialise SoS program
prog = spotsosprog;
a = msspoly('a',N);               % Define state vector
prog = prog.withIndeterminate(a); % Define independent variables

% Define system dynamics
linear_term = L_ij*a;
if (N_ijk == 0)
    f = linear_term + B_i;
else
    f = NaaSPOT(N_ijk,a) + linear_term + B_i;
end

% Define average of interest
phi = magnitude(a);

% Create V: exponents of monomials, monomials, coefficients and V
Vpow = my_monpowers(N,2*d-1);
np   = size(Vpow,1);
if abs(symmetries) == 1
    disp("BoundSys: Valid symmetry.")
    keep = prod(repmat(symmetries,np,1).^Vpow,2)==1;
    Vpow = Vpow(keep,:);
    sym = 1; % Symmetry flag 
else
    sym = 0;
end
Vm = recomp(a, Vpow, speye(size(Vpow,1)) );
[prog,Vc] = prog.newFree(size(Vpow,1));
[prog,c] = prog.newFree(1); % coefficient of energy^d
p = Vc'*Vm;
energy = a'*a;
E = c*(a'*a)^d;
V = E + p;

% f*grad(V) = f*grad(E) + f*grad(P)
fgradV = 2*d*c*energy^(d-1)*a'*(linear_term + B_i);
fgradV = fgradV + diff(p,a) * f;

% Variable to be bounded
[prog,U] = prog.newFree(1);

% SOS constraint
if epsilon == 0
    D = U - phi - fgradV;
else
    D = U - phi - fgradV - epsilon*divSPOT(sigma*sigma'*(diff(V,a))',a);
end
    
% Candidate monomials for SOS decomposition
powers = my_monpowers(N,d);
if  sym == 1
        % Exploit symmetry: 2 blocks
        np = size(powers,1);
        S1 = prod(repmat(symmetries,np,1).^powers,2)==1;
        ppp{1} = powers(S1,:);
        h1 = recomp(a, ppp{1}, speye(size(ppp{1},1)) );
        ppp{2} = powers(~S1,:);
        h2 = recomp(a, ppp{2}, speye(size(ppp{2},1)) );
        monoms = {h1; h2};
else
        h = recomp(a, powers, speye(size(powers,1)) );
        monoms = {h}; 
end  

% Minimise upper bound
switch bound
    case "U"
        % Set constraints (Upper bound)
        prog = prog.withSparseSOS(D,monoms);        
        sol = prog.minimize(U,solver,opts);
        % Compute D(a) - h^TQh
        if sym == 0
            % Non-symmetric case
            res =   D - ...
                    sol.gramMonomials{1}{1}'*double(sol.eval(sol.gramMatrices{1}{1}))*sol.gramMonomials{1}{1};
        else
            % Symmetric case
            res =   D - ...
                    sol.gramMonomials{1}{1}'*double(sol.eval(sol.gramMatrices{1}{1}))*sol.gramMonomials{1}{1} - ...
                    sol.gramMonomials{1}{2}'*double(sol.eval(sol.gramMatrices{1}{2}))*sol.gramMonomials{1}{2};
        end
            
    case "L"
        % Set constraints (Lower bound)
        prog = prog.withSparseSOS(-D,monoms);        
        sol = prog.minimize(-U,solver,opts); 
        % Compute D(a) - h^TQh
        if sym == 0
            % Non-symmetric case
            res =   -D - ...
                    sol.gramMonomials{1}{1}'*double(sol.eval(sol.gramMatrices{1}{1}))*sol.gramMonomials{1}{1};
        else
            % Symmetric case
            res =   -D - ...
                    sol.gramMonomials{1}{1}'*double(sol.eval(sol.gramMatrices{1}{1}))*sol.gramMonomials{1}{1} - ...
                    sol.gramMonomials{1}{2}'*double(sol.eval(sol.gramMatrices{1}{2}))*sol.gramMonomials{1}{2};
        end
        
    otherwise
        disp("BoundSys ERROR: First argument not valid.");
        disp("bound = 'U' for upper bound and bound = 'L' for lower bound.");
        return
        
end

U   = double(sol.eval(U)); 
% Compute residual
res = sol.eval(res);
res = norm(res.coeff,Inf);

    end
