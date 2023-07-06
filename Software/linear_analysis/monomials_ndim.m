function [Zbig]=monomials_ndim(v_vars,v_deg,n)
% This function creates a basis for polynomials of degree d in 
% n dimensions in variables 'vars'.
%
% n - dimension of basis
% v_deg - vector of degrees of polynomial entries
% v_vars - vector of variables in the expression
%
% Practically, this creates a matrix of the form
% [Z   ]
% [ Z  ]
% [  Z ]
% [   Z]
%
% The size of the resultant matrix is n*(d+m /take d) where m is the number
% of variables
%
% where the number of copies of Z is n
% version 1.0   M. Peet, Stanford. mmpeet@stanford.edu


Z=monomials(v_vars,v_deg);
Zbig=Z;
for i=2:n
    Zbig=blkdiag(Z,Zbig);
end
