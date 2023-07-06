function [prog,P]=sosposmatr(prog,n,th)
% This program declares a symbolic positive scalar semidefinite matrix P of size nxn 
% inputs: 
% prog
% n - size of matrix variable
% th - some random symbolic variable used in prog(doesn't matter which)
%
% version 1.0   M. Peet, Stanford. mmpeet@stanford.edu

Z=monomials(th,[1:n]);
[prog,VAR] = sossosvar(prog,Z,'wscoeff');
sizett=prog.var.idx{length(prog.var.idx)}-prog.var.idx{length(prog.var.idx)-1};
begin_n=prog.var.idx{length(prog.var.idx)-1};
nTT=sqrt(sizett);
nn=begin_n;
for j=1:nTT
    for i=1:nTT
        P(i,j)=sym(['coeff_',int2str(nn)]);
        nn=nn+1;
    end
end
% sizett^2=length(Z)=n=nTT
