%%%%% KSE Uncertain System %%%%%
% Example of usage of BoundUSys()
% to find a lower bound on long-time avegared energy in 
% N-diemnsional Kuramoto-Sivashinsky equation with noise.
% The system is previusly reduced to a n-dimensional uncertain
% system with n<N.

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

clear, clc;

N = 12;                      % Dimension of the N-dimensional truncated KSE
L = 2.4;                    % Length scale in KSE (Domain: [-pi*L,pi*L])
rescaling = sqrt(2*pi*L);   % Rescaling factor

% Arguments for BoundSys
magnitude = @(a,q2) (rescaling^2)*(a'*a + 2*q2)/(2*pi*L); % Magnitude to be bounded
bound = 'U';        % Looking for lower bound
d = 1;              % 2d = 4th degree auxiliary functional
n = 6;              % Dimension of the n-dimensioanl uncertain system
g = [0.02,1,0.02];  % g = [g1 g2 g3] required by UODESys 
UODESys_file = "data/KSEoutputN"+N+"n"+n+"L"+L+".mat";

%% Create the Finite-Dimensional System for KSE
mkdir("data");
f  = "data/KSEinputN" +N+"L"+L+".mat";
% Build the N-dimensional system if not done yet
if not(isfile(f))  initKSE(L,N,rescaling,f); end

%% Find an Upper Bound
[U,res,~] = BoundUSys(bound,f,{n,g,UODESys_file},magnitude,d);

disp("------------------------------");
disp("Upper Bound  : " + U);
disp("Residual norm: " + res);
disp("------------------------------");

beep;
