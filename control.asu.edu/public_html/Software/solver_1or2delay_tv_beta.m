% PD control over a distance
%cd c:\user\mmpeet\matlab\krasov_linear\linear_solver
clear all; close all; echo off; maple restart
 
%syms x1 x2 x1th x2th x1ksi x2ksi x1d1 x2d1 th ksi
syms x1 x2 x1th x2th x1ksi x2ksi x1d1 x2d1 th ksi ss hdot h

% tau=1.7177;
orderth = 2;
orderhdot=2;
otau=2;
n_dim=2;
ordernu=0;
% fact=1;
% A0=[-2 0; 0 -.9];
% A1=[-1 0;-1 -1];
 A0=[0 1; -1 -2];
 A1=[0 0;-1 1];
% tau=1.7177;
% A0=0;
% A1=-1;
 tau = h;
 taumax = 6.5;
 taumin=3
 mu=.1
fact=1;%/tau;
% data using orderth=4, ordernu=0
%mu=.01, taumin=.01, taumax=1.54
%mu=.1, taumin=.1, taumax=1.41
%mu=.3, taumin=.1, taumax=1.11
%mu=.5, taumin=.1, taumax=.81
%mu=.7, taumin=.01, taumax=.52
%mu=.75, taumin=.01, taumax=.44
%mu=.77, taumin=.01, taumax=.41
%mu=.78, taumin=.25, taumax=.4
%mu=.8, taumin=.21, taumax=.36
%mu=.82, taumin=.26, taumax=.33
%mu=.83, taumin=.26, taumax=.31
%mu=.84, taumin=.22, taumax=.27
%mu=.85, taumin=.25, taumax=.26
%mu=.86, taumin=.25, taumax=.25

%mu=.78, taumin=.01, taumax=.39
%mu=.8, taumin=.01, taumax=.34
%mu=.82, taumin=.01, taumax=.28
%mu=.83, taumin=.01, taumax=.24
%mu=.84, taumin=.01, taumax=.19
%mu=.85, taumin=.01, taumax=.14
%mu=.86, taumin=.01, taumax=.06
%mu=.865, taumin=.01, taumax=.01


% local positivity regions:
%g1=th*(th+tau);            %negative on interval [-\tau,0]
%g1=th*(th+taumax);            %negative on interval [-taumax,0]
g1=th*(th+h);            %negative on interval [-taumax,0]
g2=(hdot-mu)*(hdot+mu);
g3=(h-taumin)*(h-taumax);            %negative on interval [0,taumax]
prog = sosprogram([x1,x2,x1th,x2th,x1ksi,x2ksi,x1d1,x2d1,th,ksi,ss,hdot,h]);
% This is the single delay case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% phase 1 - variables %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we will need to declare variables:
% P, Q(th), S(th), R(th,ksi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step1 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we first declare variables P,Q1,Q2,S1,S2 and constrain positivity
disp('Creating Lyapunov Variables')
[prog,P] = sossymmatrvar(prog,sym(1),n_dim);
[prog,Q] = sosmatrvar(prog,monomials([th],0:orderth),[n_dim,n_dim]);
[prog,S] = sossymmatrvar(prog,monomials([th],0:orderth),[n_dim]);

%Now declare first spacing functions
disp('Using Spacing Functions')
[prog,U] = sossymmatrvar(prog,monomials([th,h],0:orderth),[n_dim]);
%[prog,U] = sossymmatrvar(prog,monomials([fact*th,h],0:orderth),[n_dim]);
prog = sosmatreq(prog,int(U,th,-h,0));

% now constrain P,Q,S in H1
% [P,Q1;Q1^T S1]+[U1 0;0 0] \ge 0
% [P,Q2;Q2^T S2]+[U2 0;0 0] \ge 0
if n_dim==1
    vartable1=[x1 x1th];
elseif n_dim==2
    vartable1=[x1 x2 x1th x2th];
else
    disp('oops, dimension too high')
end
    %vartable1=[x1 x1th];
%zeros22=sym(zeros(2,2));
poly1=vartable1*([P+U tau*Q;tau*Q.' tau*S])*vartable1.';

% local positivity multipliers, s1 and s2
if orderth>0
      [prog,s1] = sossosvar(prog,kron(monomials(vartable1,1),monomials([th,h],0:ceil(orderth/2))));
      [prog,s5] = sossosvar(prog,kron(monomials(vartable1,1),monomials([th,h],0:ceil(orderth/2))));
%      [prog,s1] = sossosvar(prog,kron(monomials(h,0:(ceil(orderhdot/2))),kron(monomials(vartable1,1),monomials([fact*th],0:ceil(orderth/2)-1))));
%      [prog,s5] = sossosvar(prog,kron(monomials(h,0:(ceil(orderhdot/2)-1)),kron(monomials(vartable1,1),monomials([fact*th],0:ceil(orderth/2)))));
else
   s1 = 0;
   s5 = 0;
end

eps1=.01;
% constraint P,Q,S in H1 strictly
con1 = poly1+s1*g1+s5*g3-eps1*x1^2;
prog = sosineq(prog,con1,'sparsemultipartite',{vartable1,[th,h]});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step2 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we first construct the variable R
% First declare a monomial basis for R
% disp('Using Mercer Kernels')
% Z=monomials([fact*th],0:ceil(orderth/2));
% nZ=length(Z);
% % this condition enforces global positivity of M
% %%%%%%%%%%%%%%%%%%%% hack to creat matrix P>0 %%%%%%%%%%%%%%%%%%
% % This creates a positive semidefinite matrix variable of size 4nZ*4nZ, where
% % nZ is the length of the monomial basis
% [prog,L]=sosposmatr(prog,2*nZ,th);
% 
% % Now equate the block os the matrix to pieces of Rij
% R=blkdiag(Z,Z).'*L*subs(blkdiag(Z,Z),th,ksi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step3 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we first construct the variable R
% First declare a monomial basis for R
disp('Using Semi-Seperable Kernels')
%Z=monomials([th],0:((orderth-ordernu)/2));
Z=monomials([th],0:((orderth)/2));
%Z=monomials([fact*th],0:((orderth)/2));
nZ=length(Z);
bigZ1=[];
for i=1:2*n_dim
    bigZ1=blkdiag(bigZ1,Z);
end
% This creates a positive semidefinite matrix variable of size 4nZ*4nZ, where
% nZ is the length of the monomial basis

[prog,Lsep]=sosposmatrvar(prog,2*n_dim*nZ,ordernu,ss);

LeftZ{1}=bigZ1;
RightZ{1}=subs(LeftZ{1},th,ksi);
mat_tmp=LeftZ{1}.'*Lsep*RightZ{1};
mat_tmp1=mat_tmp(1:n_dim,1:n_dim);
mat_tmp3=mat_tmp(n_dim+1:2*n_dim,1:n_dim);
mat_tmp2=mat_tmp(1:n_dim,n_dim+1:2*n_dim);
mat_tmp4=mat_tmp(n_dim+1:2*n_dim,n_dim+1:2*n_dim);
disp('integration')

% kernp1=(int(mat_tmp1,ss,-taumax,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0));
% kernp2=(int(mat_tmp1,ss,-taumax,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0));
kernp1=(int(mat_tmp1,ss,-h,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0));
kernp2=(int(mat_tmp1,ss,-h,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% phase 2 - derivatives %%%%%%%%%%%%%%%%%%%%
%%%%%%% Step1 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the derivative matrices are spelled out in the paper, they are 
disp('constructing Derivatives')
zz2=sym(zeros(n_dim,n_dim));

D= ...
[P*A0+A0.'*P+subs(Q+Q.'+S,th,0)         P*A1+subs(-Q,th,-tau)*(1-hdot)          tau*(subs(kernp1,[th ksi],[0 th])+A0.'*Q-diff(Q,th));
 (P*A1+(1-hdot)*subs(-Q,th,-tau)).'              (1-hdot)*subs(-S,th,-tau)                tau*(subs(-kernp2,[th ksi],[-tau th])*(1-hdot) + A1.'*Q);
 tau*(subs(kernp2,[ksi],[0])+A0.'*Q-diff(Q,th)).'     tau*(subs(-kernp1,[ksi],[-tau])*(1-hdot) + A1.'*Q).'     -tau*diff(S,th) ];

%  D= ...
%  [P*A0+A0.'*P+subs(Q+Q.'+S,th,0)         P*A1+subs(-Q,th,-tau)*(1-hdot)          tau*(subs(kernp1,[th ksi],[0 th])+A0.'*Q-diff(Q,th))+hdot*Q.';
%   (P*A1+(1-hdot)*subs(-Q,th,-tau)).'              (1-hdot)*subs(-S,th,-tau)                tau*(subs(-kernp2,[th ksi],[-tau th])*(1-hdot) + A1.'*Q);
%   tau*(subs(kernp2,[ksi],[0])+A0.'*Q-diff(Q,th)).'+hdot*Q.'     tau*(subs(-kernp1,[ksi],[-tau])*(1-hdot) + A1.'*Q).'     -tau*diff(S,th)+hdot*S ];


%Q11 = (A0.'*P11+P11*A0+subs(P12,th,0)+subs(P12.',th,0)+subs(P22,th,0));
%Q12 = (P11*A1+subs(P12,th,-h)*(doth-1));
%Q13 = h*(A0.'*P12-diff(P12,th)+subs(R,ksi,0));
%Q22 = (doth-1)*subs(P22,th,-h);
%Q23 = h*(A1.'*P12+(doth-1)*subs(R,ksi,-h));
%Q33 = -h*(diff(P22,th));

% D= ...
% [P*A0+A0.'*P+subs(Q+Q.'+S,th,0)         P*A1+subs(-Q,th,-tau)          tau*(subs(kernp1,[th ksi],[0 th])+A0.'*Q-diff(Q,th));
%  (P*A1+subs(-Q,th,-tau)).'              subs(-S,th,-tau)                tau*(subs(-kernp2,[th ksi],[-tau th]) + A1.'*Q);
%  tau*(subs(kernp2,[ksi],[0])+A0.'*Q-diff(Q,th)).'     tau*(subs(-kernp1,[ksi],[-tau]) + A1.'*Q).'     2*tau*(subs(kernp1-kernp2,ksi,th))-tau*diff(S,th) ];

% D= ...
% [P*A0+A0.'*P+subs(Q+Q.'+S,th,0)         P*A1+subs(-Q,th,-tau)          tau*(A0.'*Q-diff(Q,th));
%  (P*A1+subs(-Q,th,-tau)).'              subs(-S,th,-tau)                tau*(  A1.'*Q);
%  tau*(A0.'*Q-diff(Q,th)).'     tau*(A1.'*Q).'     -tau*diff(S,th) ];

kernD1=diff(kernp1,th)+diff(kernp1,ksi);
kernD2=diff(kernp2,th)+diff(kernp2,ksi);
%G=(diff(R,th)+diff(R,ksi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step2 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce positivity of D 
disp('enforcing negativity of derivative')
%Now declare first spacing functions
disp('Using Spacing Functions')
%[prog,T11] = sossymmatrvar(prog,kron(monomials(h,0:orderhdot),kron(monomials([fact*th],0:orderth),monomials([hdot],0:orderhdot))),[2*n_dim]);
[prog,T11] = sossymmatrvar(prog,monomials([fact*th,h,hdot],0:orderth),[2*n_dim]);
for i = 1:2*n_dim
   for j = i:2*n_dim
         prog = soseq(prog,int(T11(i,j),th,-h,0));
   end
end
T=blkdiag(T11,sym(zeros(n_dim,n_dim)));
% now constrain P,Q,S in H1
% [P,Q1;Q1^T S1]+[U1 0;0 0] \ge 0
% [P,Q2;Q2^T S2]+[U2 0;0 0] \ge 0
%vartable2=[x1 x2 x1d1 x2d1 x1th x2th];

if n_dim==1
    vartable2=[x1 x1d1 x1th];
elseif n_dim==2
    vartable2=[x1 x2 x1d1 x2d1 x1th x2th];;
else
    disp('oops, dimension too high')
end

%zeros22=sym(zeros(2,2));
poly2=-vartable2*(D+T)*vartable2.';
disp('creating S-variables')
% local positivity multipliers, s1 and s2
if orderth>0
%       [prog,s2] = sossosvar(prog,kron(monomials(vartable2,1),monomials([th,h,hdot],0:ceil(orderth/2))));
%       [prog,s3] = sossosvar(prog,kron(monomials(vartable2,1),monomials([th,h,hdot],0:ceil(orderth/2))));
%       [prog,s4] = sossosvar(prog,kron(monomials(vartable2,1),monomials([th,h,hdot],0:ceil(orderth/2))));
      [prog,s2] = sossosvar(prog,monomials_mp({vartable2,[th,h,hdot]},[1,ceil(orderth/2)-1]));
      [prog,s3] = sossosvar(prog,monomials_mp({vartable2,[th,h,hdot]},[1,ceil(orderth/2)-1]));
      [prog,s4] = sossosvar(prog,monomials_mp({vartable2,[th,h,hdot]},[1,ceil(orderth/2)-1]));

else
   s2 = 0;
   s3 = 0;
   s4 = 0;
end

eps2=0;%.1;
% constraint P,Q,S in H1 strictly
con2 = poly2+s2*g1+s3*g2+s4*g3-eps2*x1^2;

disp('creating multipartite constraints')
prog = sosineq(prog,con2,'sparsemultipartite',{vartable2,[th,hdot,h]});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step3 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce positivity of G
% disp('Using Mercer Kernels')
% [prog,L2]=sosposmatr(prog,2*nZ,th);
% 
% % Now equate the block os the matrix to pieces of M
% M=blkdiag(Z,Z).'*L2*subs(blkdiag(Z,Z),th,ksi);
% vartableL=[x1th, x2th];
% vartableR=[x1ksi,x2ksi];
% 
% pmin=vartableL*(M-G)*vartableR.';
% prog = soseq( prog , pmin);

disp('Using Semi-Seperable Kernels')

% % This creates a positive semidefinite matrix variable of size 4nZ*4nZ, where
% % nZ is the length of the monomial basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% New %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,L2sep]=sosposmatrvar(prog,2*n_dim*nZ,ordernu,[ss]);
% 
% % LeftZ{1}=bigZ1;
% % RightZ{1}=subs(LeftZ{1},th,ksi);
%  mat_tmp=LeftZ{1}.'*L2sep*RightZ{1};
%  mat_tmp1=mat_tmp(1:n_dim,1:n_dim);
%  mat_tmp3=mat_tmp(n_dim+1:2*n_dim,1:n_dim);
%  mat_tmp2=mat_tmp(1:n_dim,n_dim+1:2*n_dim);
%  mat_tmp4=mat_tmp(n_dim+1:2*n_dim,n_dim+1:2*n_dim);
% 
%  kernDp1=int(mat_tmp1,ss,-taumax,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0);
%  kernDp2=int(mat_tmp1,ss,-taumax,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [prog,L2sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({[ss,h,hdot]},[ceil(ordernu/2)]));
 [prog,L3sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({[ss,h,hdot]},[ceil(ordernu/2)]));
 [prog,L4sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({[ss,h,hdot]},[ceil(ordernu/2)]));
%  [prog,L2sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({ss,[h,hdot]},[ceil(ordernu/2),ceil(otau/2)]));
%  [prog,L3sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({ss,[h,hdot]},[ceil(ordernu/2),ceil(otau/2)]));
%  [prog,L4sep]=sosposmatrvar2(prog,2*n_dim*nZ,monomials_mp({ss,[h,hdot]},[ceil(ordernu/2),ceil(otau/2)]));


% LeftZ{1}=bigZ1;
% RightZ{1}=subs(LeftZ{1},th,ksi);
 mat_tmp=LeftZ{1}.'*L2sep*RightZ{1};
 mat_tmp1=mat_tmp(1:n_dim,1:n_dim);
 mat_tmp3=mat_tmp(n_dim+1:2*n_dim,1:n_dim);
 mat_tmp2=mat_tmp(1:n_dim,n_dim+1:2*n_dim);
 mat_tmp4=mat_tmp(n_dim+1:2*n_dim,n_dim+1:2*n_dim);
 kernL2_1=int(mat_tmp1,ss,-tau,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0);
 kernL2_2=int(mat_tmp1,ss,-tau,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0);
%  kernL2_1=int(mat_tmp1,ss,-tau,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0);
%  kernL2_2=int(mat_tmp1,ss,-tau,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0);

 mat_tmp=LeftZ{1}.'*L3sep*RightZ{1};
 mat_tmp1=mat_tmp(1:n_dim,1:n_dim);
 mat_tmp3=mat_tmp(n_dim+1:2*n_dim,1:n_dim);
 mat_tmp2=mat_tmp(1:n_dim,n_dim+1:2*n_dim);
 mat_tmp4=mat_tmp(n_dim+1:2*n_dim,n_dim+1:2*n_dim);
 kernL3_1=int(mat_tmp1,ss,-tau,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0);
 kernL3_2=int(mat_tmp1,ss,-tau,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0);

 mat_tmp=LeftZ{1}.'*L4sep*RightZ{1};
 mat_tmp1=mat_tmp(1:n_dim,1:n_dim);
 mat_tmp3=mat_tmp(n_dim+1:2*n_dim,1:n_dim);
 mat_tmp2=mat_tmp(1:n_dim,n_dim+1:2*n_dim);
 mat_tmp4=mat_tmp(n_dim+1:2*n_dim,n_dim+1:2*n_dim);
 kernL4_1=int(mat_tmp1,ss,-tau,th)+int(mat_tmp3,ss,th,ksi)+int(mat_tmp4,ss,ksi,0);
 kernL4_2=int(mat_tmp1,ss,-tau,ksi)+int(mat_tmp2,ss,ksi,th)+int(mat_tmp4,ss,th,0);

 % 
 vartableL=[x1th, x2th];
 vartableR=[x1ksi,x2ksi];
% 
%  pmin1=vartableL*(kernDp1-kernD1)*vartableR.';
%  pmin2=vartableL*(kernDp2-kernD2)*vartableR.';
  pmin1=vartableL*((kernL2_1-g2*kernL3_1-g3*kernL4_1)-kernD1)*vartableR.';
  pmin2=vartableL*((kernL2_2-g2*kernL3_2-g3*kernL4_2)-kernD2)*vartableR.';
disp('running equalities')
 prog = soseq( prog , pmin1);
 prog = soseq( prog , pmin2);
 
 % 
% vartableL=[x1th, x2th];
% vartableR=[x1ksi,x2ksi];
% 
% pmin1=vartableL*(kernDp1-kernD1)*vartableR.';
% pmin2=vartableL*(kernDp2-kernD2)*vartableR.';
% disp('running equalities')
% prog = sosmatreq( prog , kernDp1-kernD1);
% prog = sosmatreq( prog , kernDp1-kernD1);


disp('Computing Solution')

pars.alg=2;
pars.stepdif=1;
pars.eps=10^(-18);
pars.maxiter=100;
pars.cg.maxiter=600;
pars.cg.qprec=1;
pars.cg.stagtol=1e-22;
pars.cg.restol=5e-5;
%pars.free=0;

prog = sossolve(prog);
%prog = sossolve(prog,pars);

% for i=1:length(Lsep)
%     for j=1:length(Lsep)
%         Lsepn(i,j)=sosgetsol(prog,Lsep(i,j));
%     end
% end
% 
% for i=1:length(kernp1)
%     for j=1:length(kernp1)
%         kernp1n(i,j)=sosgetsol(prog,kernp1(i,j));
%     end
% end

%prog = sossolve(prog);
% for i=1:length(L2)
%     for j=1:length(L2)
%         L2n(i,j)=sosgetsol(prog,L2(i,j));
%     end
% end
% for i=1:length(L)
%     for j=1:length(L)
%         Ln(i,j)=sosgetsol(prog,L(i,j));
%     end
% end

%prog = sossolve(prog);

