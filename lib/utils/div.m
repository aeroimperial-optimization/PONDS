function D = div(v,x)
% Compute divergence of v with respect to x when working with YALMIP
% sdpavar

% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

N = size(v,1);
D = 0;
for i = 1:N
    D = D + jacobian(v(i),x(i));
end