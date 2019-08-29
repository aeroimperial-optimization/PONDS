function [P, Pc, Pm] = symmpoly(x, DEGREE, SYMMETRIES)

% [P, Pc, Pm] = SYMMPOLY(x, DEGREE, SYMMETRIES) returns an N-variate tunable 
%       polynomial P(x) of degree DEGREE, which is invariant under the M 
%       sign symmetries specified by the rows of the M-by-N matrix. That is,
%
%       P( SYMMETRIES(i,:) .* x ) = P(x)
%
%       The tunable coefficients of P(x) are returned in the vector Pc,
%       while the corresponding monomials are listed in Pm.
%
%       Inputs:
%       - x: an N-by-1 (or 1-by-N) sdpvar
%       - SYMMETRIES: an M-by-N matrix with entries +1 or -1
%
%       Example. Build a degree-4 polynomial P(x,y) invariant under the
%       transformations (x,y) -> (-x,y) and (x,y) -> (x,-y).
%
%       >> sdpvar x y
%       >> SYMMETRIES = [-1, 1; 1, -1];
%       >> [P, C] = symmpoly([x,y], 4, SYMMETRIES);

% Get all monomial powers
N = length(x);
powers = monpowers(N, DEGREE);

% Find invariant monomials
keep = prod(SYMMETRIES(1,:).^powers, 2) == 1;       
for i = 2:size(SYMMETRIES, 1)
   keep = keep & prod(SYMMETRIES(i,:).^powers, 2) == 1;
end
powers  = powers(keep, :);

% Build P
Pc = sdpvar(size(powers, 1), 1);
Pm = recovermonoms(powers, x);
P  = Pc' * Pm;