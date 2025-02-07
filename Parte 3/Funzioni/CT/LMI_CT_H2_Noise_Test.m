function [K,rho,feas]=LMI_CT_H2_Noise_Test(A,B,C,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system,
% including white noise and H2 norm constraint.

Btot=[];
for i=1:N
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);
mtot=sum(m);

% White noise covariance matrix
R_sqrt= eye(18); % Adjust size according to system dimensions

% H2 constraint matrices
C_H2 = [eye(36); zeros(18,36)];
D_H2 = [zeros(36,18); R_sqrt];

% White noise input
W = eye(18); % White noise

% Optimization variables
S = sdpvar(mtot+ntot,mtot+ntot);

if ContStruc==ones(N,N)
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
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

% LMI Constraints
LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'+eye(ntot)+ W*W' <=-1e-2*eye(ntot), Y>=1e-2*eye(ntot)];

% H2 Constraint
H2_constr = [S C_H2*Y + D_H2*L; L'*D_H2' + Y*C_H2' Y];
n = size(H2_constr, 1);
H2_constrain = (H2_constr >= 1e-2 * eye(n));

% Objective function
objective_function = trace(S);

% Solve optimization problem
constrains = [LMIconstr, H2_constrain];
options = sdpsettings('solver','sedumi');
J = optimize(constrains, objective_function, options);
feas = J.problem;
L = double(L);
Y = double(Y);

K = L/Y;
rho = max(real(eig(A+Btot*K)));
end
