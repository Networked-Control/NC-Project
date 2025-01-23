function [K,rho,feas]=LMI_CT_Effort(A,B,C,N,ContStruc,alpha_L,alpha_Y)
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
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);   % ntot = 36
assignin('base', 'ntot', ntot);
mtot=sum(m);      % mtot = 18
assignin('base', 'mtot', mtot);

k_L =sdpvar;
k_Y = sdpvar;

if ContStruc==ones(N,N)
    % Centralized design
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        Y=blkdiag(Y,sdpvar(n(i)));
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

    LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'<=-1e-2*eye(ntot)]+[Y>=1e-2*eye(ntot)];

    M = [sqrt(k_L)*eye(ntot) L';L eye(mtot)];  % 54x54
    assignin('base', 'M', M);

    constr_1 = [M>=1e-2*eye(ntot+mtot)];

    M_2 = [k_Y*eye(ntot) eye(ntot);eye(ntot) Y]; % 72x72
    assignin('base', 'M_2', M_2);

    constr_2 = [M_2>=1e-2*eye(ntot*2)];

%Optimization

constraints = [constr_1,constr_2,LMIconstr];
objective_function = alpha_L*k_L+alpha_Y*k_Y;
options=sdpsettings('solver','sedumi');

J = optimize(constraints,objective_function,options);

L= double(L);
Y= double(Y);

feas=J.problem;

K=L/Y;
rho=max(real(eig(A+Btot*K)));
