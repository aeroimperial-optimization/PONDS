%%%%% Lorenz Attractor %%%%%
% Compute dissipation(rho) in Lorenz attractor using BoundSys()
% and plot the results

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

clear, clc;
mkdir("data");

% Value of the parameters of Lorenz attractor
rho = 0:5:50;
beta  = 8/3; 
sigma = 10;
rescaling = 20; % Rescaling factor

% Arguments for BoundSys
bound = 'U';                    % Looking for upper bound
d = 4;                          % 2d = 8th degree auxiliary functional
epsilon = 0; sigma_noise = 0;   % No noise
verbose = 1;                    % Enable verbosity for SDP solver
symmetries = 0;                 % No simetries enable
SOSPsolver = 'spotless';        % SOSP solver
SDPsolver = "mosek";            % SDP solver

% Magnitude to be bounded rescaled to the original value
magnitude = @(a) (rescaling^2)*(sigma*(a(1))^2 + a(2)^2 + beta*a(3)^2);

U = zeros(1,length(rho));
res = zeros(1,length(rho));
for i = 1:length(rho)
    % Create the 3-simensional system for Lorenz attractor
    f  = "data/LORENZinput"+"beta"+beta+"sigma"+sigma+"rho"+rho(i)+".mat";
    % Build the initial system if not done yet
    if not(isfile(f))  initLorenz(beta,sigma,rho(i),rescaling,f); end
    % Call BoundSys
    [U(i), res(i),~]  = BoundSys(bound,f,magnitude,d,epsilon,sigma_noise,verbose,symmetries,SOSPsolver,SDPsolver);
end

%% Plot U and res
subplot(2,1,1)
plot(rho,U,"bo--");
grid on, xlabel("$\rho$","Interpreter","latex"), ylabel("Dissipation","Interpreter","latex");
subplot(2,1,2)
plot(rho,res,"bo--");
grid on, xlabel("$\rho$","Interpreter","latex"), ylabel("Residual","Interpreter","latex");


