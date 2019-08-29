function v = Naa(N_ijk,a)
% Compute N_ijk(a)a
%
% N_ijk     Third order tensor of size (N,N,N)
% a         vector of size N
%
% Written by Mario Lino Valencia (September 2019)
% Imperial College London - Department of Aeronautics

N = length(a);
aa = a*a';
v = double2sdpvar(zeros(N,1));
for i = 1:N
    N_i  = N_ijk(i,:,:);
    if (N_i == 0) continue; end 
    v(i) = N_i(:)'*aa(:);
end
