function [prog,P] = sossymmatrvar(prog,Z,sizeP);
%
% version 1.0   A. Papachristodoulou
%               M. Peet

for i = 1:sizeP*(sizeP+1)/2
   [prog,matrPi(i)] = sospolyvar(prog,Z,'wscoeff');
end

dumvec(1) = 1;
for i = 1:sizeP-1
   dumvec(i+1) = dumvec(i)+(sizeP-i)+1;
end
dumvec(end+1) = dumvec(end)+1;

dumQ = diag([-1:-1:-sizeP]);
for i = 1:sizeP
   dumQ = dumQ + diag(dumvec(i):dumvec(i+1)-1,i-1)+ diag(dumvec(i):dumvec(i+1)-1,-i+1);
end

for i = 1:sizeP;
   for j = 1:sizeP;
      P(i,j) = matrPi(dumQ(i,j));
   end
end
