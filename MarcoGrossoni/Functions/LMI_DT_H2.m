function [K,rho,feas]=LMI_DT_H2(F,G,H,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - A: system matrix.
% - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)
Btot=[];

for i=1:N
    Gdec{i}=G(:,i);
    Hdec{i}=(H(2*(i-1)+1:2*i,:))';
end

Gtot=[];
for i=1:N
    m(i)=size(Gdec{i},2);
    n(i)=size(Hdec{i},2);
    Gtot=[Gtot,Gdec{i}];
end
ntot=size(F,1);
mtot=sum(m);

R_sqrt= [1 0 0;
        0 1 0;
        0 0 1];

H_H2 = [eye(6);zeros(3,6)];
D_H2 = [zeros(6,3);R_sqrt];


yalmip clear

S = sdpvar(mtot+ntot,mtot+ntot);

if ContStruc==ones(N,N)
    % Centralized design
    P=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Dentralized/distributed design
    P=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        P=blkdiag(P,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end
end

LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F' Gtot*L;
        L'*Gtot' P]>=1e-1*eye(ntot*2)];
H2_constr = [S H_H2*P + D_H2*L; L'*D_H2' + P*H_H2' P] >= 1e-2*(eye(15));

objective_function = trace(S);

constrains = [LMIconstr,H2_constr];
options=sdpsettings('solver','sdpt3');
J=optimize(constrains,objective_function,options);
feas=J.problem;
L=double(L);
P=double(P);

K=L/P;
rho=max(real(eig(F+Gtot*K)));


