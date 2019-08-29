function v = Nab(N_ijk,a,b)
% Compute N_ijk(a)a
%
% N_ijk     Third order tensor of size (N,Na,Nb)
% a         vector of dimension Na
% b         vector of dimension Nb
%
% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

N  = size(N_ijk,1);
ab = a*b';
v = double2sdpvar(zeros(N,1));
for i = 1:N
    N_i  = N_ijk(i,:,:);
    if (N_i == 0) continue; end 
    v(i) = N_i(:)'*ab(:);
end
