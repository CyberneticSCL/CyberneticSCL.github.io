function [prog,Q]=sosposmatrvar(prog,n,d,vars)
% This program declares a symbolic positive scalar semidefinite matrix P of size nxn 
% inputs: 
% prog
% n - dimension of matrix
% d - degree of polynomial entries
% vars - variables in the expression
%
% creates a matrix variable of the form
% [Z   ]^T   [Z   ]
% [ Z  ]   Q [ Z  ]
% [  Z ]     [  Z ]
% [   Z]     [   Z], Q>0
%
% version 1.0   M. Peet, Stanford. mmpeet@stanford.edu


Z=monomials(vars,0:ceil(d/2));
Zbig=Z;
for i=2:n
    Zbig=blkdiag(Z,Zbig);
end
nZbig=size(Zbig,1);

[prog,P]=sosposmatr(prog,nZbig,vars(1));
Q=Zbig.'*P*Zbig;

