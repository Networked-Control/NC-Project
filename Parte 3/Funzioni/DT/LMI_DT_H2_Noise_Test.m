function [K, rho, feas] = LMI_DT_H2_Noise_Test(F, G, H, N, ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system,
% considering white noise constraints.
% 
% Inputs:
% - F: system matrix.
% - G: input matrices.
% - H: output matrices.
% - N: number of subsystems.
% - ContStruc: NxN matrix specifying the control information structure.
%
% Outputs:
% - K: structured control gain.
% - rho: spectral abscissa of (F + G*K).
% - feas: feasibility of the LMI problem (0 if feasible).

Gtot = [];
for i = 1:N
    m(i) = size(G{i}, 2);
    n(i) = size(H{i}, 1);
    Gtot = [Gtot, G{i}];
end
ntot = size(F, 1);
mtot = sum(m);

R_sqrt = eye(18);
H_H2 = [eye(36); zeros(18, 36)];
D_H2 = [zeros(36, 18); R_sqrt];

yalmip clear
S = sdpvar(mtot+ntot, mtot+ntot);

if ContStruc == ones(N, N)
    P = sdpvar(ntot);
    L = sdpvar(mtot, ntot);
else
    P = [];
    L = sdpvar(mtot, ntot);
    minc = 0;
    for i = 1:N
        P = blkdiag(P, sdpvar(n(i)));
        ninc = 0;
        for j = 1:N
            if ContStruc(i, j) == 0
                L(minc+1:minc+m(i), ninc+1:ninc+n(j)) = zeros(m(i), n(j));
            end
            ninc = ninc + n(j);
        end
        minc = minc + m(i);
    end
end

% White noise covariance matrix
Q_w = eye(ntot); % Identity matrix for white noise

% LMI Constraints including noise effect
LMIconstr = [[P - F*P*F' - F*L'*Gtot' - Gtot*L*F' + Q_w, Gtot*L;
              L'*Gtot', P] >= 1e-1*eye(ntot*2)];
          
H2_constr = [S, H_H2*P + D_H2*L;
             L'*D_H2' + P*H_H2', P];
         
[n, ~] = size(H2_constr);
H2_constrain = (H2_constr >= 1e-2 * eye(n));

objective_function = trace(S);

constrains = [LMIconstr, H2_constr];
options = sdpsettings('solver', 'sdpt3');
J = optimize(constrains, objective_function, options);
feas = J.problem;
L = double(L);
P = double(P);

K = L / P;
rho = max(real(eig(F + Gtot*K)));
end
