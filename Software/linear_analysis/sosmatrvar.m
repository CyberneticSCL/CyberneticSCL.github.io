function [prog,P] = sosmatrvar(prog,Z,sizePQ);

sizeP = sizePQ(1);
sizeQ = sizePQ(2);

for i = 1:sizeP*sizeQ
   [prog,matrPi(i)] = sospolyvar(prog,Z,'wscoeff');
end

for i = 1:sizeP;
   for j = 1:sizeQ;
      P(i,j) = matrPi((i-1)*sizeQ+j);
   end
end
