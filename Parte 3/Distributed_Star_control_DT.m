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
% Discrete time 
Autovalori_DT = eig(F);
isStable_2 = true;
for i = 1:length(Autovalori_DT)
    if abs(real(Autovalori_DT(i))) > 1+1e-4
        isStable_2 = false;
    end
end
if isStable_2
    disp('System is Stable in DT');
else
    disp('System is Unstable in DT');
end

moduli = abs(Autovalori_DT);
spectral_radius = max(moduli);  % Spectral Radius
disp(['Spectral Radius: ', num2str(spectral_radius)]);

%% Control Structures
alpha = 1;  % Must be positive, the negative sign is already considered in the LMI computation
%rho_DT = exp(alpha*Ts);
rho_DT = 0.88;
center = -0.6; % Must be positive, the negative sign is already considered in the LMI computation
radius = 0.45; % center and radius are computed for Circle LMIs
angle = 45; % Sector LMIs
alpha_L = 0.1; % Effort LMIs
alpha_Y = 10; % Effort LMIs

% Distributed LMI Performance
ContStruc_Distr_star=eye(N);
for i=2:N
    ContStruc_Distr_star(1,i)=1;
    ContStruc_Distr_star(i,1)=1;
end
[star_fm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Distr_star,3);

% Discrete Time
[K_star_DT,rho_star_DT,feas_star_DT]=LMI_DT_Stability(F,Gd,Hd,N,ContStruc_Distr_star); % LMI for stability
[K_star_DT_perf,rho_star_DT_perf,feas_star_DT_perf]=LMI_DT_Performance(F,Gd,Hd,N,ContStruc_Distr_star,rho_DT); % LMI for performance
[K_star_DT_circle,rho_star_CT_circle,feas_star_CT_circle]=LMI_DT_Circle_Area(A,Bd,Cd,N,ContStruc_Distr_star,center,radius) % LMI for circle delimited area
[K_star_DT_effort,rho_star_DT_effort,feas_star_DT_effort]=LMI_DT_Effort(F,Gd,Hd,N,ContStruc_Distr_star,alpha_L,alpha_Y);
[K_star_DT_H2,rho_star_DT_H2,feas_star_DT_H2]=LMI_DT_H2(F,Gd,Hd,N,ContStruc_Distr_star);

%% Display
 disp('Results (Discrete-time):')
 disp(['-  Distributed_Star_DT: Feasibility=',num2str(feas_star_DT),', rho=',num2str(rho_star_DT),', FM=',num2str(star_fm_DT),'.'])
 disp(['-  Distributed_Star_DT_Perf: Feasibility=',num2str(feas_star_DT_perf),', rho=',num2str(rho_star_DT_perf),', FM=',num2str(star_fm_DT),'.'])
 disp(['-  Distributed_Star_DT_Circle: Feasibility=',num2str(feas_star_CT_circle),', rho=',num2str(rho_star_CT_circle),', FM=',num2str(star_fm_DT),'.'])
 disp(['-  Distributed_Star_DT_Effort: Feasibility=',num2str(feas_star_DT_effort),', rho=',num2str(rho_star_DT_effort),', FM=',num2str(star_fm_DT),'.'])
 disp(['-  Distributed_Star_DT_H2: Feasibility=',num2str(feas_star_DT_H2),', rho=',num2str(rho_star_DT_H2),', FM=',num2str(star_fm_DT),'.'])

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
% DT Simulation

Q = 0.01 * eye(n); % Varianza del rumore bianco
a = 0.1; % Definisce l'ampiezza del rumore

k = 0;
for k=1:Tfinal/Ts
    w = a * (2 * rand(n,1) - 1); % rumore uniforme

    % state variable
    x_star_DT(:,k)=((F+G*K_star_DT)^k)*x0;
    x_star_DT_perf(:,k)=((F+G*K_star_DT_perf)^k)*x0;
    x_star_DT_circle(:,k)=((F+G*K_star_DT_circle)^k)*x0;
    x_star_DT_effort(:,k)=((F+G*K_star_DT_effort)^k)*x0;
    x_star_DT_H2(:,k)=((F+G*K_star_DT_H2)^k)*x0 + w;

    % control variable
    u_star_DT(:,k) = K_star_DT * x_star_DT(:,k);
    u_star_DT_perf(:,k) = K_star_DT_perf * x_star_DT_perf(:,k);
    u_star_DT_circle(:,k) = K_star_DT_circle * x_star_DT_circle(:,k);
    u_star_DT_effort(:,k) = K_star_DT_effort * x_star_DT_effort(:,k);
    u_star_DT_H2(:,k) = K_star_DT_H2 * x_star_DT_H2(:,k);
end

%% Calcolo autovalori
eig_DT = eig(F+G*K_star_DT_circle)

%% Creazione della figura
figure;
hold on;
grid on;
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

% Disegna il cerchio unitario
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5); % Cerchio unitario

% Disegna il cerchio desiderato
plot(-center + radius * cos(theta), radius * sin(theta), 'r-', 'LineWidth', 1.5); 

% Disegna gli autovalori
plot(real(eig_DT), imag(eig_DT), 'bx', 'MarkerSize', 10, 'LineWidth', 2);

% Disegna gli assi
plot([-1.2, 1.2], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-1.2, 1.2], 'k', 'LineWidth', 1);

% Titolo e legenda
title('Autovalori e Regione di Collocazione');
legend({'Cerchio unitario', 'Regione desiderata', 'Autovalori'}, 'Location', 'Best');
xlabel('Re');
ylabel('Im');

hold off;


%% Plotting px1 coordinate for every LMI used
figure
plot([Ts:Ts:Tfinal],x_star_DT(1,:),[Ts:Ts:Tfinal],x_star_DT_perf(1,:),[Ts:Ts:Tfinal],x_star_DT_circle(1,:),[Ts:Ts:Tfinal],x_star_DT_effort(1,:),[Ts:Ts:Tfinal],x_star_DT_H2(1,:))
title('DT controller Position in X')
grid on
legend('DT Stability', 'DT Performance', 'DT Circle Area', 'DT Effort', 'DT H2')
xlabel('Time (k)')
ylabel('Position (X)')

%% Plotting py1 coordinate for every LMI used
figure
plot([Ts:Ts:Tfinal],x_star_DT(3,:),[Ts:Ts:Tfinal],x_star_DT_perf(3,:),[Ts:Ts:Tfinal],x_star_DT_circle(3,:),[Ts:Ts:Tfinal],x_star_DT_effort(3,:),[Ts:Ts:Tfinal],x_star_DT_H2(3,:))
title('DT controller Position in Y')
grid on
legend('DT Stability', 'DT Performance', 'DT Circle Area', 'DT Effort', 'DT H2')
xlabel('Time (k)')
ylabel('Position (Y)')

%% Plotting ux control effort for every LMI used
figure
plot([Ts:Ts:Tfinal],u_star_DT(1,:),[Ts:Ts:Tfinal],u_star_DT_perf(1,:),[Ts:Ts:Tfinal],u_star_DT_circle(1,:),[Ts:Ts:Tfinal],u_star_DT_effort(1,:),[Ts:Ts:Tfinal],u_star_DT_H2(1,:))
title('DT control variable in X')
grid on
legend('DT Stability', 'DT Performance', 'DT Circle Area', 'DT Effort', 'DT H2')
xlabel('Time (k)')
ylabel('Control U(x)')

%% Plotting uy control effort for every LMI used
figure
plot([Ts:Ts:Tfinal],u_star_DT(2,:),[Ts:Ts:Tfinal],u_star_DT_perf(2,:),[Ts:Ts:Tfinal],u_star_DT_circle(2,:),[Ts:Ts:Tfinal],u_star_DT_effort(2,:),[Ts:Ts:Tfinal],u_star_DT_H2(2,:))
title('DT control variable in Y')
grid on
legend('DT Stability', 'DT Performance', 'DT Circle Area', 'DT Effort', 'DT H2')
xlabel('Time (k)')
ylabel('Control U(y)')