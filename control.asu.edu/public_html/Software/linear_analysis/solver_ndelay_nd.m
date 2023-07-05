%cd c:\user\mmpeet\matlab\krasov_linear\example_2delay
clear; close all; echo off; maple restart
% This program determines stablity of a linear differential equation with a
% a single delay, where \dot{x}(t)=A0x(t)+A{1}x(t-tau(1))+...+A{K}x(t-tau(K))  
% where A0, A{i}, and tau(i) are user inputs. 
%
% Inputs: A{i} - these can be arbitrary square matrices of arbitrary 
%         dimension. However,  the higher the higher the dimension of A{i},
%         the more time the program will take to run
%         
%         tau(i) - These can be an arbitrary sequence of positive increasing 
%         numbers.
%
%         orderth - This input controls the accuracy of the results. For
%         most problems, orderth=2 should be sufficient to obtain a
%         reasonable degree of accuracy. Note: orderth should be an even
%         integer.
% 
% Requirements: In order to operate, this program requires a working
%               version of SOStools. There are some known compatability 
%               issues with SOStools and Matlab version 7+ due to errors in 
%               the implementation of Maple v8. In addition, it
%               is highly recommended that the user have a compiled version
%               of cdd. Finally the following package of subprograms are
%               required which allows SOStools to handle matrix objects
%               directly:
%               SOStools Matrix Addon package:
%                   sosposmatr.m
%                   sosmatrvar.m
%                   sossymmatrvar.m
%                   sosposmatrvar.m
%
% version 1.0   M. Peet, Stanford. mmpeet@stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%
% Enter system dynamics in terms of A0, A{1}, ..., A{K}
% A0=[0 1; 0 -.1];
% A{1}=[0; .1]*[-3.75 -11.5];
%====================
% Gu 2delay example 1
% A0=[-2 0; 0 -.9];
% A{1}=[-1 0;-1 -1]*.05;
% A{2}=[-1 0;-1 -1]*.95;
% tau(2) = 8.5
% tau(1) = tau(2)/2;
%====================
% Gu 2delay example 2
% A0=[0 1; -1 .1];
% A{1}=[0 0;-1 0];
% A{2}=[0 0;1 0];
% tau(2) = 1.3722
% tau(1) = tau(2)/2;

% Enter values of the delay tau(1), ..., tau(K)


% tau(1) = 1.166;
%tau(2) = 1.34;
%
% Enter degree of accuracy - must be an even integer
orderth = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fact=1/tau;
fact=1/5;           % conditioning multiplier
eps1=.01;           % strictness of lyapunov positivity
eps2=.01;           % strictness of derivative negativity

% control inputs to SeDuMi 
pars.alg=2;
pars.stepdif=1;
pars.eps=10^(-10);
pars.maxiter=100;
pars.cg.maxiter=200;
pars.cg.qprec=1;
pars.cg.stagtol=1e-22;
pars.cg.restol=5e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal processing:
syms th ksi

n_dim=length(A0);
n_delay=length(tau);
tauK=tau(n_delay);

% declare lots of symbolic variables
vartablex=[];
vartablexth=[];
vartablexksi=[];
for j=1:n_delay
    eval(['vartabled',int2str(j),'=[];'])
end   
for i=1:n_dim
    eval(['syms x',int2str(i),' x',int2str(i),'th x',int2str(i),'ksi']);
    eval(['vartablex=[vartablex, x',int2str(i),'];'])
    eval(['vartablexth=[vartablexth, x',int2str(i),'th];'])
    eval(['vartablexksi=[vartablexksi, x',int2str(i),'ksi];'])
    for j=1:n_delay
        eval(['syms x',int2str(i),'d',int2str(j)]);
        eval(['vartabled',int2str(j),'=[vartabled',int2str(j),', x',int2str(i),'d',int2str(j),'];'])
    end        
end
mastervartable=[th,ksi,vartablex,vartablexth,vartablexksi];
for j=1:n_delay
    eval(['mastervartable=[mastervartable, vartabled',int2str(j),'];'])
end 

prog = sosprogram(mastervartable);

% local positivity regions:
g{1}=th*(th+tau(1));            %negative on interval [-\tau,0]
for i=2:n_delay
    g{i}=(th+tau(i))*(th+tau(i-1));%negative on interval [-tau_i,-\tau_i-1]
end

% This is the multiple delay case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% phase 1 - variables %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we will need to declare variables:
% P, Q(th), S(th), R(th,ksi)
%
disp('Creating Lyapunov Variables')
[prog,P] = sossymmatrvar(prog,sym(1),[n_dim]);
for i=1:n_delay
    eval(['[prog,Q',int2str(i),'] = sosmatrvar(prog,monomials([fact*th],0:orderth),[n_dim,n_dim]);'])
    eval(['[prog,S',int2str(i),'] = sossymmatrvar(prog,monomials([fact*th],0:orderth),[n_dim]);'])
end

%Now declare first spacing functions
disp('Using Spacing Functions')
for i=1:n_delay
    eval(['[prog,U',int2str(i),'] = sossymmatrvar(prog,monomials([fact*th],0:orderth),[n_dim]);'])
end

U_constr=int(vartablex*(U1)*vartablex.',th,-tau(1),0);
for i=2:n_delay
    eval(['U_constr=U_constr+int(vartablex*(U',int2str(i),')*vartablex.'',th,-tau(',int2str(i),'),-tau(',int2str(i-1),'));'])
end

% Implementation 1
prog = soseq(prog,U_constr);


% now constrain P,Q,S in H1
vartable1=[vartablex vartablexth];

for i=1:n_delay
    eval(['poly',int2str(i),'=vartable1*([P+U',int2str(i),' tau(n_delay)*Q',int2str(i),';tau(n_delay)*Q',int2str(i),'.'' tau(n_delay)*S',int2str(i),'])*vartable1.'';'])
end
%poly2=vartable1*([P+U2 tau*Q2;tau*Q2.' tau*S2])*vartable1.';

% local positivity multipliers, s1 and s2
zzn=sym(zeros(n_dim,n_dim));

if orderth>0
    for i=1:n_delay
        eval(['[prog,bigs',int2str(i),'] = sosposmatrvar(prog,2*n_dim,orderth-2,[th]);'])
    end
else
    for i=1:n_delay
        eval(['bigs',int2str(i),'=blkdiag(zzn,zzn);'])
    end
end


% constraint P,Q,S in H1 strictly

extras=vartablex*vartablex.';

for i=1:n_delay
    eval(['tempcon = -eps1*extras+poly',int2str(i),'+vartable1*(bigs',int2str(i),'*g{',int2str(i),'})*vartable1.'';'])
    prog = sosineq(prog,tempcon,'sparsemultipartite',{vartable1,[th]});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step2 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we first construct the variable R
% First declare a monomial basis for R
disp('Using Mercer Kernels')
%Z=monomials([fact*th],0:ceil(orderth/2));
Z=monomials([fact*th],0:0*orderth);
nZ=length(Z);
bigZ1=[];
for i=1:n_dim
    bigZ1=blkdiag(bigZ1,Z);
end
% This creates a positive semidefinite matrix variable of size 4nZ*4nZ, where
% nZ is the length of the monomial basis
[prog,L]=sosposmatr(prog,n_dim*n_delay*nZ,th);

% Now equate the block os the matrix to pieces of Rij
delta(1) = tau(1);
for i=2:n_delay 
    delta(i)=tau(i)-tau(i-1);
end

LeftZ{1}=subs(bigZ1,th,tauK/delta(1)*th+0*tauK/delta(1));
RightZ{1}=subs(LeftZ{1},th,ksi);
for i=2:n_delay
    delta(i)=tau(i)-tau(i-1);
    LeftZ{i}=subs(bigZ1,th,tauK/delta(i)*th+tau(i-1)*tauK/delta(i));
    RightZ{i}=subs(LeftZ{i},th,ksi);
end

for i=1:n_delay
    for j=1:n_delay
        R{i,j}=LeftZ{i}.'*L((i-1)*nZ*n_dim+1:i*n_dim*nZ,(j-1)*nZ*n_dim+1:j*n_dim*nZ)*RightZ{j};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% phase 2 - derivatives %%%%%%%%%%%%%%%%%%%%
%%%%%%% Step1 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the derivative matrices are spelled out in the paper, they are 
disp('constructing Derivatives')
zzn=sym(zeros(n_dim,n_dim));
D11=sym(zeros(n_dim*(n_delay+1),n_dim*(n_delay+1)));
D11(1:n_dim,1:n_dim)=P*A0+A0.'*P+subs(Q1+Q1.'+S1,th,0);
eval(['D11(1:n_dim,n_delay*n_dim+1:(n_delay+1)*n_dim)=P*A{n_delay}-subs(Q',int2str(n_delay),',th,-tau(n_delay));'])
eval(['D11(n_delay*n_dim+1:(n_delay+1)*n_dim,1:n_dim)=(P*A{n_delay}-subs(Q',int2str(n_delay),',th,-tau(n_delay))).'';'])
eval(['D11(n_delay*n_dim+1:(n_delay+1)*n_dim,n_delay*n_dim+1:(n_delay+1)*n_dim)=-subs(S',int2str(n_delay),',th,-tau(n_delay)).'';'])
for i=1:(n_delay-1)
    eval(['D11(n_dim*i+1:n_dim*(i+1),n_dim*i+1:n_dim*(i+1))=subs(S',int2str(i+1),'-S',int2str(i),',th,-tau(i));'])
    eval(['D11(1:n_dim,(i*n_dim+1):((i+1)*n_dim))=P*A{i}+subs(Q',int2str(i+1),'-Q',int2str(i),',th,-tau(i));'])
    eval(['D11((i*n_dim+1):((i+1)*n_dim),1:n_dim)=(P*A{i}+subs(Q',int2str(i+1),'-Q',int2str(i),',th,-tau(i))).'';'])
end
for i=1:n_delay
    eval(['Qi=Q',int2str(i),';'])
    % define D_12{i}
    D_temp=[];
    for j=1:(n_delay-1)
        D_temp=[D_temp; subs(subs(R{j+1,i}-R{j,i},th,-tau(j)),ksi,th) + A{j}.'*Qi];
    end
    D12i=[subs(subs(R{1,i},th,0),ksi,th)+A0.'*Qi-diff(Qi,th);
        D_temp;
        A{n_delay}.'*Qi-subs(subs(R{n_delay,i},th,-tau(n_delay)),ksi,th) ];

    eval(['D22i=-diff(S',int2str(i),',th);'])

    Dt{i}=[D11 tauK*D12i;tauK*D12i.' tauK*D22i];
end

for i=1:n_delay
    for j=1:n_delay
        G{i,j}=(diff(R{i,j},th)+diff(R{i,j},ksi));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step2 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce negativity of D 
disp('enforcing negativity of derivative')
%Now declare first spacing functions
disp('Using Spacing Functions')

vartablexxd=[vartablex];
for i=1:n_delay
    eval(['vartablexxd=[vartablexxd vartabled',int2str(i),'];'])
end

for i=1:n_delay
    eval(['[prog,T',int2str(i),'] = sossymmatrvar(prog,monomials([fact*th],0:orderth),[(n_delay+1)*n_dim]);'])
end
T_constr=int(vartablexxd*(T1)*vartablexxd.',th,-tau(1),0);
for i=2:n_delay
    eval(['T_constr=T_constr+int(vartablexxd*(T',int2str(i),')*vartablexxd.'',th,-tau(',int2str(i),'),-tau(',int2str(i-1),'));'])
end
% Implementation 1
prog = soseq(prog,T_constr);

% now constrain -D_i+T_i in H1
vartable2=[vartablexxd vartablexth];
for i=1:n_delay
   eval(['Tfull{',int2str(i),'}=blkdiag(T',int2str(i),',zzn);'])
   newpoly{i}=vartable2*(Tfull{i}-Dt{i})*vartable2.';
end
if orderth>0
    for i=1:n_delay
        eval(['[prog,newbigs',int2str(i),'] = sosposmatrvar(prog,(n_delay+2)*n_dim,orderth-2,[th]);'])
    end
else
    for i=1:n_delay
        eval(['newbigs',int2str(i),'=sym(zeros((n_delay+2)*n_dim,(n_delay+2)*n_dim));'])
    end
end
for i=1:n_delay
%    eval(['tempcon = newpoly',int2str(i),'+vartable2*(newbigs',int2str(i),'*g{',int2str(i),'}-eps2*sym(blkdiag(eye(n_dim),zeros((n_delay+1)*n_dim,(n_delay+1)*n_dim))))*vartable2.'';'])
    eval(['tempcon = -eps2*extras+newpoly{i}+vartable2*(newbigs',int2str(i),'*g{',int2str(i),'})*vartable2.'';'])
    prog = sosineq(prog,tempcon,'sparsemultipartite',{vartable2,[th]});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step3 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce positivity of G
disp('Using Mercer Kernels')

% Now equate the block os the matrix to pieces of M
vartableL=vartablexth;%=[x1th, x2th];
vartableR=vartablexksi;%=[x1ksi,x2ksi];

% Now equate the block os the matrix to pieces of M
[prog,L2]=sosposmatr(prog,n_delay*n_dim*nZ,th);
for i=1:n_delay
    for j=1:n_delay
        Mij=LeftZ{i}.'*L2((i-1)*nZ*n_dim+1:i*n_dim*nZ,(j-1)*nZ*n_dim+1:j*n_dim*nZ)*RightZ{j};
        pmin_temp=vartableL*(Mij-G{i,j})*vartableR.';
        prog = soseq( prog , pmin_temp);
    end
end

disp('Computing Solution')
prog = sossolve(prog,pars);



