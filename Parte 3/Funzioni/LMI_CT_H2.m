function [K,rho,feas]=LMI_CT_H2(A,B,C,N,ContStruc)
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
ntot=size(A,1);
mtot=sum(m);

R_sqrt= [1 0 0;
        0 1 0;
        0 0 1];

% C_H2 = [eye(6);zeros(3,6)];
% disp(C_H2);
% D_H2 = [zeros(6,3);R_sqrt];
% disp(D_H2)

n = 54; % Numero di righe
C_H2 = rand(n, 36); % Matrice random n x 36
D_H2 = rand(n, 18); % Matrice random n x 18


S = sdpvar(mtot+ntot,mtot+ntot);

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

LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'+eye(ntot)<=-1e-2*eye(ntot)]+[Y>=1e-2*eye(ntot)];

%disp(size(S));            % 54 x 54
%disp(L);                  % 18 x 36
%disp(size(Y));            % 36 x 36 

H2_constr = [S C_H2*Y + D_H2*L; L'*D_H2' + Y*C_H2' Y];

% Determina la dimensione della matrice H2_constr
[n, ~] = size(H2_constr); % n è il numero di righe

% Costruisci il vincolo con la matrice identità della stessa dimensione
H2_constrain = (H2_constr >= 1e-2 * eye(n));

objective_function = trace(S);

constrains = [LMIconstr,H2_constrain];
options=sdpsettings('solver','sdpt3');
J=optimize(constrains,objective_function,options);
feas=J.problem;
L=double(L);
Y=double(Y);

K=L/Y;
rho=max(real(eig(A+Btot*K)));


