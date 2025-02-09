close all 
clear all
clc

addpath('Funzioni/CT');
addpath('Funzioni/DT');
addpath('Funzioni');

%% MATRIXES DEFINITION

%Array of masses
% Marcello Farina, 22/10/2020

N=9;
m=[3 2 3 2 3 2 4 1 2];
k=[0 0.1 0 0.2 0 0 0 0 0
   0.1 0 0.2 0 0.1 0 0 0 0
   0 0.2 0 0 0 0.3 0 0 0
   0.2 0 0 0 0.2 0 0.5 0 0
   0 0.1 0 0.2 0 0.7 0 0.1 0
   0 0 0.3 0 0.7 0 0 0 0.4
   0 0 0 0.5 0 0 0 0.9 0
   0 0 0 0 0.1 0 0.9 0 0.1
   0 0 0 0 0 0.4 0 0.1 0];

h=[0 0.13 0 0.13 0 0 0 0 0
   0.13 0 0.2 0 0.14 0 0 0 0
   0 0.2 0 0 0 0.3 0 0 0
   0.13 0 0 0 0.2 0 0.1 0 0
   0 0.14 0 0.2 0 0.2 0 0.1 0
   0 0 0.3 0 0.2 0 0 0 0.4
   0 0 0 0.1 0 0 0 0.1 0
   0 0 0 0 0.1 0 0.1 0 0.4
   0 0 0 0 0 0.4 0 0.4 0];

A=zeros(N*4,N*4);
n = N*2;
B=[];
for i=1:N
    for j=1:N
        if i==j
            Ai=[0 1;-sum(k(i,1:N))/m(i) -sum(h(i,1:N))/m(i)];
        else
            Ai=[0 0;k(i,j)/m(i) h(i,j)/m(i)];
        end
        Ac{i,j}=blkdiag(Ai,Ai);
        A((i-1)*4+1:4*i,(j-1)*4+1:4*j)=Ac{i,j};
    end
    B=blkdiag(B,blkdiag([0 1/m(i)]',[0 1/m(i)]'));
end
C=eye(4*N);

%% Discretization
% Time costant analysis
Ts= 0.1; %Shannon Theorem on signal sampling (fs=2*fc at least) 
[F,G,H,W,Ts]=ssdata(c2d(ss(A,B,C,[]),Ts));

%% Decomposition

for i=1:2:(N*2-1)
    Bd{ceil(i/2)}=B(:,i:i+1);
    Gd{ceil(i/2)}=G(:,i:i+1);
end

for i=1:4:N*4
    Cd{ceil(i/4)}=C(i:i+3,:);
    Hd{ceil(i/4)}=H(i:i+3,:);
end

%% Stability
% Continuous time
Autovalori_CT = eig(A);
isStable = true;
for i = 1:length(Autovalori_CT)
    if real(Autovalori_CT(i)) > 0
        isStable = false;
    end
end
if isStable
    disp('System is Stable in CT');
else
    disp('System is Unstable in CT');
end

real_parts_CT = real(Autovalori_CT);

spectral_abscissa = max(real_parts_CT);   % Spectral Abscissa
disp(['Spectral Abscissa: ', num2str(spectral_abscissa)]);

%% Control Structures
alpha = 1;
%rho_DT = exp(alpha*Ts);
rho_DT = 0.88;
center = 20; % Must be positive, the negative sign is already considered in the LMI computation
radius = 1; % center and radius are computed for Circle LMIs
angle = 45; % Sector LMIs
alpha_L = 10; % Effort LMIs
alpha_Y = 0; % Effort LMIs


% Distributed LMI Performance
ContStruc_Distr_star=eye(N);
for i=2:N
    ContStruc_Distr_star(1,i)=1;
    ContStruc_Distr_star(i,1)=1;
end
[star_fm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Distr_star,3);

% Continuous Time 
[K_star_CT,rho_star_CT,feas_star_CT]=LMI_CT_Stability(A,Bd,Cd,N,ContStruc_Distr_star); % LMI for stability
[K_star_CT_perf,rho_star_CT_perf,feas_star_CT_perf]=LMI_CT_Performance(A,Bd,Cd,N,ContStruc_Distr_star,alpha); % LMI for performance
[K_star_CT_sector,rho_star_CT_sector,feas_star_CT_sector]=LMI_CT_Sector(A,Bd,Cd,N,ContStruc_Distr_star,angle) % LMI for sector delimited area
[K_star_CT_effort,rho_star_CT_effort,feas_star_CT_effort]=LMI_CT_Effort(A,Bd,Cd,N,ContStruc_Distr_star,alpha_L,alpha_Y) % LMI for sector delimited area
[K_star_CT_H2,rho_star_CT_H2,feas_star_CT_H2]=LMI_CT_H2(A,Bd,Cd,N,ContStruc_Distr_star) % LMI for H2

%% Display

 disp('Results (Continuous-time):')
 disp(['-  Distributed_Star_CT: Feasibility=',num2str(feas_star_CT),', rho=',num2str(rho_star_CT),', FM=',num2str(star_fm),'.'])
 disp(['-  Distributed_Star_CT_Perf: Feasibility=',num2str(feas_star_CT_perf),', rho=',num2str(rho_star_CT_perf),', FM=',num2str(star_fm),'.'])
 disp(['-  Distributed_Star_CT_Sector: Feasibility=',num2str(feas_star_CT_sector),', rho=',num2str(rho_star_CT_sector),', FM=',num2str(star_fm),'.'])
 disp(['-  Distributed_Star_CT_Effort: Feasibility=',num2str(feas_star_CT_effort),', rho=',num2str(rho_star_CT_effort),', FM=',num2str(star_fm),'.'])
 disp(['-  Distributed_Star_CT_H2: Feasibility=',num2str(feas_star_CT_H2),', rho=',num2str(rho_star_CT_H2),', FM=',num2str(star_fm),'.'])

%% Plots
Gtot=[];
Htot=[];
Btot=[];
Ctot=[];
for i=1:N
    Btot=[B,Bd{i}];
    Ctot=[C
        Cd{i}];
    Gtot=[G,Gd{i}];
    Htot=[H
        Hd{i}];
end

% simulation data
Tfinal=10;
T=0:0.01:Tfinal;
% Random initial condition 
min_x0 = -10;
max_x0 = 10;
x0 = min_x0 + (max_x0-min_x0) .* rand(36,1);  % Generate a random number between [min_x0, max_x0]

k = 0;
Q = eye(size(G,2));
x_star_CT_H2(:,1) = x0;
a = 0.1; % Definisce l'ampiezza del rumore

% CT Simulation 
for t=T
    k=k+1;
     w = a * (2 * rand(2*n,1) - 1); % rumore uniforme

    % state computation
    x_star_free(:,k)=expm(A*t)*x0; % No control
    x_star_CT(:,k)=expm((A+B*K_star_CT)*t)*x0;
    x_star_CT_perf(:,k)=expm((A+B*K_star_CT_perf)*t)*x0;
    x_star_CT_sector(:,k)=expm((A+B*K_star_CT_sector)*t)*x0;
    x_star_CT_effort(:,k)=expm((A+B*K_star_CT_effort)*t)*x0;
    x_star_CT_H2(:,k)=expm((A+B*K_star_CT_H2)*t)*x_star_CT_H2(:,1) + w;

    % control variable
    u_star_CT(:,k) = K_star_CT * x_star_CT(:,k);
    u_star_CT_perf(:,k) = K_star_CT_perf * x_star_CT_perf(:,k);
    u_star_CT_sector(:,k) = K_star_CT_sector * x_star_CT_sector(:,k);
    u_star_CT_effort(:,k) = K_star_CT_effort * x_star_CT_effort(:,k);
    u_star_CT_H2(:,k) = K_star_CT_H2 * x_star_CT_H2(:,k);
end

% Continuous Time figure
% Primo grafico: posizione lungo X
figure
plot(T, x_star_free(1,:),T, x_star_CT(1,:), T, x_star_CT_perf(1,:), T, x_star_CT_sector(1,:),T,x_star_CT_effort(1,:),T,x_star_CT_H2(1,:)) % Position of the first Mass along x direction
title('CT controller Position in X')  
grid on
legend('No control','CT Stability', 'CT Performance', 'CT Sector', 'CT Effort', 'CT H2') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (X)') % Etichetta dell'asse y

% Secondo grafico: posizione lungo Y
figure
plot(T, x_star_free(3,:),T, x_star_CT(3,:), T, x_star_CT_perf(3,:), T, x_star_CT_sector(3,:),T, x_star_CT_effort(3,:),T, x_star_CT_H2(3,:)) % Position of the first Mass along y direction
title('CT controllers Position in Y')
grid on
legend('No control','CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (Y)') % Etichetta dell'asse y

% Terzo grafico: variabile di controllo U lungo X
figure
plot(T, u_star_CT(1,:),T, u_star_CT_perf(1,:),T, u_star_CT_sector(1,:),T, u_star_CT_effort(1,:),T, u_star_CT_H2(1,:)) % variabile di controllo lungo x
title('CT control variable in X')
grid on
legend('CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Control action (U(x))') % Etichetta dell'asse x

figure
plot(T, u_star_CT(2,:),T, u_star_CT_perf(2,:),T, u_star_CT_sector(2,:),T, u_star_CT_effort(2,:),T, u_star_CT_H2(2,:)) % variabile di controllo lungo y
title('CT control variable in Y')
grid on
legend('CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse y
ylabel('Control action (U(y))') % Etichetta dell'asse y

%Eigen = eig(A-B*K_c_CT_effort)