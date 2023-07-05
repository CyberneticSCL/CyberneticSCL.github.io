function b = int_p(a,x,x1,x2)
%
% DESCRIPTION 
%   A blatent copy of diff.m to
%   Integrate a polynomial.
%   
% INPUTS 
%   A: polynomial 
%   X: integrate with respect to X.
%   x1: lower limit of integration
%   x2: upper limit of integration
%
% OUTPUTS  
%   B: polynomial
%  
% SYNTAX 
%   B = int(A,X);
%     Integrate the polynomial, A, with respect to X, from x1 to x2.  A 
%     should be a polynomial and X should be a polynomial variable or a 
%     string. x1 and x2 should be scalars or scalar polynomials.  
%     Integration is done element-by-element if A is a matrix.

% 6/20/2010: MMP  Initial Coding   -   based on diff.m by PJS

% Error Checking
if nargin~=4 
  error('Error in calling diff');
end

if isa(x,'polynomial')
  x = x.varname{1};
elseif ~ischar(x)
  error('X must be a polynomial variable or a string');
end

% Get polynomial info about a
a = polynomial(a);
a = combine(a);
adeg = a.degmat;
avar = a.varname;
nta = size(adeg,1);
[nra,nca]=size(a);
acoef = reshape(a.coefficient,nta,nra*nca);

% Get format x1 and x2
x1 = polynomial(x1);
x1 = combine(x1);
x2 = polynomial(x2);
x2 = combine(x2);

% Find variable we are integrating with respect to.
varnumb = 0;
for i1 = 1:length(avar)
  if strcmp(avar{i1},x)
    varnumb = i1;
  end
end

% Integrate
  for i1 = 1:nta 
    acoef(i1,:) = acoef(i1,:)/(adeg(i1,varnumb)+1);
  end
  adeg(:,varnumb) = adeg(:,varnumb)+1;
  a_int = combine(polynomial(acoef,adeg,avar,[nra nca]));
  b = combine(subs(a_int,x,x2)-subs(a_int,x,x1));
% end



