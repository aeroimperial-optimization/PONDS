%%%%% KSE Finite-Dimensional System %%%%%
% Example of usage of BoundSys()
% to find a lower bound on long-time avegared energy in 
% N-diemnsional Kuramoto-Sivashinsky equation with noise.

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

clear, clc;

N = 6;                      % Dimension of the N-dimensional truncated KSE
L = 1.2;                    % Length scale in KSE (Domain: [-pi*L,pi*L])
rescaling = sqrt(2*pi*L);   % Rescaling factor

% Arguments for BoundSys
magnitude = @(a) (rescaling^2)*(a'*a)/(2*pi*L); % Magnitude to be bounded
bound = 'L';                % Looking for lower bound
d = 3;                      % 2d = 10th degree auxiliary functional
epsilon = 1e-3/rescaling;   % Rescaled noise intensity
sigma   = eye(N);
verbose = 1;                % Enable verbosity for SDP solver
symmetries = (-1).^(1:N);   % Simetries to enable block diagonalisation in SDP solver 
SOSPsolver = 'spotless';    % Using SPOTLess as SOSP solver
SDPsolver = 'mosek';        % Using Mosek as SDP solver

%% Create the Finite-Dimensional System for KSE
mkdir("data");
f  = "data/KSEinputN" +N+"L"+L+".mat";
% Build the N-dimensional system if not done yet
if not(isfile(f))  initKSE(L,N,rescaling,f); end

%% Find a Lower Bound
[L,res,~] = BoundSys(bound,f,magnitude,d,epsilon,sigma,verbose,symmetries,SOSPsolver,SDPsolver);

disp("------------------------------");
disp("Lower Bound  : " + L);
disp("Residual norm: " + res);
disp("------------------------------");

beep;
