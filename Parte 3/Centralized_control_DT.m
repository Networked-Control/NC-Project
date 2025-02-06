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
center = -0.5; % Must be positive, the negative sign is already considered in the LMI computation
radius = 0.1; % center and radius are computed for Circle LMIs
angle = 45; % Sector LMIs
alpha_L = 10^100; % Effort LMIs
alpha_Y = 0; % Effort LMIs

% Centralized LMI Performance
ContStruc_Centr = ones(N,N);
[cfm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Centr,3);
[cfm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Centr,3);

 % Discrete Time
 [K_c_DT,rho_c_DT,feas_c_DT]=LMI_DT_Stability(F,Gd,Hd,N,ContStruc_Centr); % LMI for stability
 [K_c_DT_perf,rho_c_DT_perf,feas_c_DT_perf]=LMI_DT_Performance(F,Gd,Hd,N,ContStruc_Centr,rho_DT); % LMI for performance
 [K_c_DT_circle,rho_c_DT_circle,feas_c_DT_circle]=LMI_DT_Circle_Area(F,Gd,Hd,N,ContStruc_Centr,center,radius); % LMI for performance


%% Display
 disp('Results (Discrete-time):')
 disp(['-  Centralized: Feasibility=',num2str(feas_c_DT),', rho=',num2str(rho_c_DT),', FM=',num2str(cfm_DT),'.'])
 disp(['-  Centralized_Perf: Feasibility=',num2str(feas_c_DT_perf),', rho=',num2str(rho_c_DT_perf),', FM=',num2str(cfm_DT),'.'])
 disp(['-  Centralized_Circle: Feasibility=',num2str(feas_c_DT_circle),', rho=',num2str(rho_c_DT_circle),', FM=',num2str(cfm_DT),'.'])

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
min_x0 = 0;
max_x0 = 1;
random_number = min_x0 + (max_x0-min_x0) .* rand(1,1);  % Generate a random number between [min_x0, max_x0]
x0 = repmat(random_number,36,1);

% DT Simulation

k = 0;
for k=1:Tfinal/Ts
     x_c_DT(:,k)=((F+G*K_c_DT)^k)*x0;
     x_c_DT_perf(:,k)=((F+G*K_c_DT_perf)^k)*x0;
     x_c_DT_circle(:,k)=((F+G*K_c_DT_circle)^k)*x0;
end

Eigen = eig(F+G*K_c_DT_circle)

figure
% Plotting px1 coordinate for every LMI used
plot([Ts:Ts:Tfinal],x_c_DT(1,:),[Ts:Ts:Tfinal],x_c_DT_perf(1,:),[Ts:Ts:Tfinal],x_c_DT_circle(1,:))
title('DT controllers graphs')
grid on
legend('DT Stability', 'DT Performance', 'DT Circle Area')
xlabel('Time (k)')
ylabel('Position (X)')
