function [prog,Pth]=sosposmatrvar2(prog,n,Zth,th)
% This program declares a positive semidefinite matrix P(vartable) using vartable Zth of size nxn 
% inputs: 
% prog
% n - size of matrix variable
% Zth - vector of monomials contained in the decomposition Z' M Z
% th - some random symbolic variable used in prog(doesn't matter which)

nZth=length(Zth);
if nZth>0
%create positive matrix variable P of size n*nZth x n*nZth
[prog,P]=sosposmatr(prog,n*nZth,th);

% create basis for factorization Pth=DZ'*P*DZ
DZ=sym(zeros(0));
for i=1:n
    DZ=blkdiag(DZ,Zth);
end
Pth=DZ.'*P*DZ;
else 
    Pth=sym(zeros(n,n));
end