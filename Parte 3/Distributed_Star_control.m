close all 
clear all
clc

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
alpha = 2;  % Must be positive, the negative sign is already considered in the LMI computation
rho_DT = exp(alpha*Ts);
center = 20; % Must be positive, the negative sign is already considered in the LMI computation
radius = 1;

% Distributed LMI Performance
ContStruc_Distr_star=eye(N);
for i=2:N
    ContStruc_Distr_star(1,i)=1;
    ContStruc_Distr_star(i,1)=1;
end
[star_fm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Distr_star,3);
[star_fm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Distr_star,3);

% Continuous Time 
[K_star_CT,rho_star_CT,feas_star_CT]=LMI_CT_DeDicont(A,Bd,Cd,N,ContStruc_Distr_star); % LMI for stability
[K_star_CT_perf,rho_star_CT_perf,feas_star_CT_perf]=LMI_CT_DeDicont_perf(A,Bd,Cd,N,ContStruc_Distr_star,alpha); % LMI for performance
[K_star_CT_circle,rho_star_CT_circle,feas_star_CT_circle]=LMI_Circle_Area_CT(A,Bd,Cd,N,ContStruc_Distr_star,center,radius) % LMI for circle delimited area

% Discrete Time
[K_star_DT,rho_star_DT,feas_star_DT]=LMI_DT_DeDicont(F,Gd,Hd,N,ContStruc_Distr_star); % LMI for stability
[K_star_DT_perf,rho_star_DT_perf,feas_star_DT_perf]=LMI_DT_DeDicont_perf(F,Gd,Hd,N,ContStruc_Distr_star,rho_DT); % LMI for performance

%% Display

disp('Results (Continuous-time):')
disp(['-  Distributed (star): Feasibility=',num2str(feas_star_CT),', rho=',num2str(rho_star_CT),', FM=',num2str(star_fm),'.'])
disp(['-  Distributed (star)_Perf: Feasibility=',num2str(feas_star_CT_perf),', rho=',num2str(rho_star_CT_perf),', FM=',num2str(star_fm),'.'])
disp(['-  Distributed (star)_Circle: Feasibility=',num2str(feas_star_CT_circle),', rho=',num2str(rho_star_CT_circle),', FM=',num2str(star_fm),'.'])

disp('Results (Discrete-time):')
disp(['-  Distributed (star): Feasibility=',num2str(feas_star_DT),', rho=',num2str(rho_star_DT),', FM=',num2str(star_fm_DT),'.'])
disp(['-  Distributed (star)_Perf: Feasibility=',num2str(feas_star_DT_perf),', rho=',num2str(rho_star_DT_perf),', FM=',num2str(star_fm_DT),'.'])

%% Plots
% Gtot=[];
% Htot=[];
% Btot=[];
% Ctot=[];
% for i=1:N
%     Btot=[B,Bd{i}];
%     Ctot=[C
%         Cd{i}];
%     Gtot=[G,Gd{i}];
%     Htot=[H
%         Hd{i}];
% end
% 
% % simulation data
% Tfinal=10;
% T=0:0.01:Tfinal;
% % Random initial condition 
% min_x0 = 1;
% max_x0 = 3;
% random_number = min_x0 + (max_x0-min_x0) .* rand(1,1);  % Generate a random number between [min_x0, max_x0]
% x0 = repmat(random_number,36,1);
% k = 0;
% % CT Simulation 
% for t=T
%     k=k+1;
%     x_c_CT(:,k)=expm((A+B*K_De_CT)*t)*x0;
%     x_c_CT_perf(:,k)=expm((A+B*K_CT_perf)*t)*x0;
%     x_c_CT_circle(:,k)=expm((A+B*K_CT_circle)*t)*x0;
% end
% 
% % Continuous Time figure
% figure
% % Plotting px1 coordinate for every LMI used
% %plot(T,x_c_CT(1,:),T,x_c_CT_perf(1,:),T,x_c_CT_circle(1,:))
% plot(T,x_c_CT(1,:),T,x_c_CT_perf(1,:))
% title('CT controllers graphs')
% grid on



%% PLOTS TO BE ADJUSTED
% Gtot=[];
% Htot=[];
% Btot=[]; 
% Ctot=[];
% for i=1:N
%     Btot=[B,Bd{i}];
%     Ctot=[C
%         Cd{i}];
%     Gtot=[G,Gd{i}];
%     Htot=[H
%         Hd{i}];
% end
% 
% % simulation data
% Tfinal=10;
% T=0:0.01:Tfinal;
% min_x0 = 1;
% max_x0 = 3;
% random_number = min_x0 + (max_x0-min_x0) .* rand(1,1);  % Generate a random number between [min_x0, max_x0]
% x0 = repmat(random_number,36,1);
% k=0;
% for t=T
%     k=k+1;
%     x_c(:,k)=expm((A+B*K_c)*t)*x0;
% end
% for k=1:Tfinal/Ts
%     x_c_DT(:,k)=((F+G*K_De_DT)^k)*x0;
% end
% 
% figure
% for i=1:N
%     subplot(N,2,2*(i-1)+1)
%     hold on
%     grid on
%     title(['\px_{',num2str(i),'}'])
%     plot(T,[x_c((i-1)*4+1,:)],'k')
%     axis([0 T(end) min(x0)-10 max(x0)+10])
% 
%     hold on
%     grid on
%     title(['\py_{',num2str(i),'}'])
%     plot(T,[x_c((i-1)*4+3,:)],'k')
%     axis([0 T(end) min(x0)-10 max(x0)+10])
% 
% 
%     subplot(N,2,2*i)
%     hold on
%     grid on
%     title(['\px_{',num2str(i),'}'])
%     plot(Ts:Ts:Tfinal,[x_c_DT((i-1)*4+1,:)],'k.-')
%     axis([0 T(end) min(x0)-10 max(x0)+10])
% 
%     hold on
%     grid on
%     title(['\py_{',num2str(i),'}'])
%     plot(Ts:Ts:Tfinal,[x_c_DT((i-1)*4+3,:)],'k.-')
%     axis([0 T(end) min(x0)-10 max(x0)+10])
% end
% legend('Centralized')