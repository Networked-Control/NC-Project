close all 
clear all
clc

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
% Centralized
ContStruc_Centr = ones(N,N);
[cfm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Centr,3);
[cfm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Centr,3);

[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bd,Cd,N,ContStruc_Centr);
[K_c_DT,rho_c_DT,feas_c_DT]=LMI_DT_DeDicont(F,Gd,Hd,N,ContStruc_Centr);

%% Decentralized
ContStruc_Dec = diag(ones(N,1));
[Dfm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Dec,3);
[Dfm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Dec,3);

[K_De,rho_De,feas_De]=LMI_CT_DeDicont(A,Bd,Cd,N,ContStruc_Dec);
[K_De_DT,rho_De_DT,feas_De_DT]=LMI_DT_DeDicont(F,Gd,Hd,N,ContStruc_Dec);

%% Distributed (String)
ContStruc_Distr_string=eye(N);
for i=1:N-1
    ContStruc_Distr_string(i,i+1)=1;
    ContStruc_Distr_string(i+1,i)=1;
end
[string_fm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Distr_string,3);
[string_fm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Distr_string,3);

[K_string,rho_string,feas_string]=LMI_CT_DeDicont(A,Bd,Cd,N,ContStruc_Distr_string);
[K_string_DT,rho_string_DT,feas_string_DT]=LMI_DT_DeDicont(F,Gd,Hd,N,ContStruc_Distr_string);


%% Distributed (Star bi)
ContStruc_Distr_star=eye(N);
for i=2:N
    ContStruc_Distr_star(1,i)=1;
    ContStruc_Distr_star(i,1)=1;
end
[star_fm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Distr_star,3);
[star_fm_DT]=di_fixed_modes(F,Gd,Hd,N,ContStruc_Distr_star,3);

[K_star,rho_star,feas_star]=LMI_CT_DeDicont(A,Bd,Cd,N,ContStruc_Distr_star);
[K_star_DT,rho_star_DT,feas_star_DT]=LMI_DT_DeDicont(F,Gd,Hd,N,ContStruc_Distr_star);

%% Display

disp('Results (Continuous-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c),', rho=',num2str(rho_c),', FM=',num2str(cfm),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_De),', rho=',num2str(rho_De),', FM=',num2str(Dfm),'.'])
disp(['-  Distributed (string): Feasibility=',num2str(feas_string),', rho=',num2str(rho_string),', FM=',num2str(string_fm),'.'])
disp(['-  Distributed (star): Feasibility=',num2str(feas_star),', rho=',num2str(rho_star),', FM=',num2str(star_fm),'.'])

disp('Results (Discrete-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c_DT),', rho=',num2str(rho_c_DT),', FM=',num2str(cfm_DT),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_De_DT),', rho=',num2str(rho_De_DT),', FM=',num2str(Dfm_DT),'.'])
disp(['-  Distributed (string): Feasibility=',num2str(feas_string_DT),', rho=',num2str(rho_string_DT),', FM=',num2str(string_fm_DT),'.'])
disp(['-  Distributed (star): Feasibility=',num2str(feas_star_DT),', rho=',num2str(rho_star_DT),', FM=',num2str(star_fm_DT),'.'])
