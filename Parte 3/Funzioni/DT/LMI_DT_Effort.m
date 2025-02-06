function [K,rho,feas,k_L,k_Y]=LMI_DT_effort(F,G,H,N,ContStruc,alpha_L,alpha_Y)
% Computes, using LMIs, the distributed "state feedback" control law for the discrete-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - F: system matrix.
% - G: input matrices (i.e., G{1},..., G{N} are the input matrices of the decomposed system, one for each channel).
% - H: output matrices  (i.e., H{1},..., H{N} are the output matrices of the decomposed system, one for each channel, where [Hdec{1}',...,
% Hdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral radius of matrix (F+G*K) - note that [H{1}',...,
% H{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)
Gtot=[];
for i=1:N
    m(i)=size(G{i},2);
    n(i)=size(H{i},2);
    Gtot=[Gtot,G{i}];
end
ntot=size(F,1);
mtot=sum(m);

yalmip clear

k_L =sdpvar;
k_L_sqrt = sqrt(k_L);
k_Y = sdpvar;



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
                disp(['i=', num2str(i), ', j=', num2str(j)]);
                disp(['minc=', num2str(minc), ', m(i)=', num2str(m(i))]);
                disp(['ninc=', num2str(ninc), ', n(j)=', num2str(n(j))]);
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end
end

%Control effort constrains

 
    LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F' Gtot*L;
        L'*Gtot' P]>=1e-1*eye(ntot*2)];

    constr_1 = [[k_L_sqrt*eye(ntot) L';L eye(mtot)]>=1e-1*eye(ntot+mtot)];

    constr_2 = [[k_Y*eye(ntot) eye(ntot);eye(ntot) P]>=1e-2*eye(ntot*2)];

%Optimization

constraints = [LMIconstr,constr_1,constr_2];
objective_function = alpha_L*k_L+alpha_Y*k_Y;
options=sdpsettings('solver','sedumi');

% Display constraints
disp('Constraints:');
disp(sdisplay(constraints));

J=optimize(constraints,objective_function,options);

feas=J.problem;
L=double(L);
P=double(P);

K=L/P;
rho=max(abs(eig(F+Gtot*K)));
