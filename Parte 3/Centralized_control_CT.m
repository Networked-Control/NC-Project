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
alpha = 1;  % Must be positive, the negative sign is already considered in the LMI computation
%rho_DT = exp(alpha*Ts);
rho_DT = 0.88;
angle = pi/4; % Sector LMIs
alpha_L = 0.1; % Effort LMIs
alpha_Y = 10; % Effort LMIs

% Centralized LMI Performance
ContStruc_Centr = ones(N,N);
[cfm]=di_fixed_modes(A,Bd,Cd,N,ContStruc_Centr,3);

% Continuous Time 
%% Stability
 [K_c_CT,rho_c_CT,feas_c_CT]=LMI_CT_Stability(A,Bd,Cd,N,ContStruc_Centr); % LMI for stability
%% Performance
 [K_c_CT_perf,rho_c_CT_perf,feas_c_CT_perf]=LMI_CT_Performance(A,Bd,Cd,N,ContStruc_Centr,alpha); % LMI for performance
%% Sector
 [K_c_CT_sector,rho_c_CT_sector,feas_c_CT_sector]=LMI_CT_Sector(A,Bd,Cd,N,ContStruc_Centr,angle) % LMI for sector delimited area
%% Effort 
 [K_c_CT_effort,rho_c_CT_effort,feas_c_CT_effort]=LMI_CT_Effort(A,Bd,Cd,N,ContStruc_Centr,alpha_L,alpha_Y) % LMI for sector delimited area
%% H2 
 [K_c_CT_H2,rho_c_CT_H2,feas_c_CT_H2]=LMI_CT_H2_Noise_Test(A,Bd,Cd,N,ContStruc_Centr) % LMI for H2
 %% Mixed
 [K_c_CT_Mixed,rho_c_CT_Mixed,feas_c_CT_Mixed]=LMI_CT_Mixed(A,Bd,Cd,N,ContStruc_Centr,alpha,angle,alpha_L,alpha_Y) % LMI for H2

 %% Display

 disp('Results (Continuous-time):')
 disp(['-  Centralized: Feasibility=',num2str(feas_c_CT),', rho=',num2str(rho_c_CT),', FM=',num2str(cfm),'.'])
 disp(['-  Centralized_Perf: Feasibility=',num2str(feas_c_CT_perf),', rho=',num2str(rho_c_CT_perf),', FM=',num2str(cfm),'.'])
 disp(['-  Centralized_Sector: Feasibility=',num2str(feas_c_CT_sector),', rho=',num2str(rho_c_CT_sector),', FM=',num2str(cfm),'.'])
 disp(['-  Centralized_Effort: Feasibility=',num2str(feas_c_CT_effort),', rho=',num2str(rho_c_CT_effort),', FM=',num2str(cfm),'.'])
 disp(['-  Centralized_H2: Feasibility=',num2str(feas_c_CT_H2),', rho=',num2str(rho_c_CT_H2),', FM=',num2str(cfm),'.'])
 disp(['-  Centralized_Mixed: Feasibility=',num2str(feas_c_CT_Mixed),', rho=',num2str(rho_c_CT_Mixed),', FM=',num2str(cfm),'.'])

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

%% Calcolo autovalori Stability
eig_CT = eig(A+B*K_c_CT);
spectral_abscissa_stab = max(eig_CT); 

%% Plot Autovalori Stability
figure;
hold on;
grid on;
axis equal;

% Impostiamo limiti degli assi in base agli autovalori
y_limit = max(abs(imag(eig_CT))) + 0.2;
x_limit = y_limit;

% Limiti degli assi (estesi di un po' rispetto agli autovalori)
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]);

% Disegnare solo gli autovalori
h1 = plot(real(eig_CT), imag(eig_CT), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Autovalori in blu

% Disegnare gli assi
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1); % Asse X
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1); % Asse Y

% Titolo e etichette degli assi
title('LMI Stability');
xlabel('Re');
ylabel('Im');

% Legenda per gli autovalori
legend(h1, {'Eigenvalues'}, 'Location', 'Best');

hold off;


%% Calcolo autovalori Performance
eig_CT_perf = eig(A+B*K_c_CT_perf);
spectral_abscissa_perf = max(eig_CT_perf); 

%% Plot Autovalori Performance
% Creazione della figura
figure;
hold on;
grid on;
axis equal;

% Impostiamo limiti degli assi simmetrici
x_limit = max(abs(real(eig_CT_perf))) + 0.2;
y_limit = x_limit; % Rendi il plot quadrato
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]);

% Disegnare la retta verticale (perpendicolare all'asse x) che passa per alpha
h_line = plot(-alpha * ones(1, 100), linspace(-y_limit, y_limit, 100), 'r-', 'LineWidth', 1.5); % Reetta rossa perpendicolare

% Disegnare gli autovalori
h1 = plot(real(eig_CT_perf), imag(eig_CT_perf), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Autovalori in blu

% Disegnare gli assi
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1);

% Titolo e legenda, includendo solo la linea rossa e gli autovalori
title('LMI Performance');
legend([h_line, h1], {'alpha', 'Eigenvalues'}, 'Location', 'Best');
xlabel('Re');
ylabel('Im');

hold off;

%% Calcolo autovalori Sector
eig_CT_sector = eig(A+B*K_c_CT_sector);
spectral_abscissa_sector = max(eig_CT_sector); 

%% Plot Autovalori Sector

% Find the maximum absolute real value for symmetric axes
x_limit = max(abs(real(eig_CT_sector))) + 0.2;
y_limit = x_limit; % Make the plot square

% Create the figure
figure;
hold on;
grid on;
axis equal;

% Set symmetric axis limits
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]); 

% Define the sector lines
x_sector = linspace(-x_limit, 0, 100); % Generate points for the lines
y_sector1 = tan(angle) * x_sector; % First line (+angle w.r.t. x-axis)
y_sector2 = -tan(angle) * x_sector; % Second line (-angle w.r.t. x-axis)

% Plot the sector with the two boundary lines
h1 = plot(x_sector, y_sector1, 'r-', 'LineWidth', 1.5); % Positive real part
h2 = plot(x_sector, y_sector2, 'r-', 'LineWidth', 1.5); % Negative real part

% Plot the eigenvalues in continuous time
h3 = plot(real(eig_CT_sector), imag(eig_CT_sector), 'bx', 'MarkerSize', 10, 'LineWidth', 2);

% Draw the axes with symmetric length
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1);

% Title and legend with the sector angle
title('LMI Sector');
legend([h1, h3], {sprintf('Sector (angle = %d°)', rad2deg(angle)), 'Eigenvalues'}, 'Location', 'Best');
xlabel('Re');
ylabel('Im');

hold off;

%% Calcolo autovalori Effort
eig_CT_effort = eig(A+B*K_c_CT_effort);
spectral_abscissa_effort = max(eig_CT_effort);   % Spectral Abscissa

%% Plot Autovalori Effort
% Definizione dei valori di alpha_L e alpha_Y
figure;
hold on;
grid on;
axis equal;

% Impostiamo limiti degli assi in base agli autovalori
y_limit = max(abs(imag(eig_CT_effort))) + 0.2;
x_limit = y_limit;

% Limiti degli assi (estesi di un po' rispetto agli autovalori)
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]);

% Disegnare solo gli autovalori
h1 = plot(real(eig_CT_effort), imag(eig_CT_effort), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Autovalori in blu

% Disegnare gli assi
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1); % Asse X
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1); % Asse Y

% Titolo e etichette degli assi
title('LMI Effort');
xlabel('Re');
ylabel('Im');

% Creare oggetti fittizi per alpha_L e alpha_Y, così non avranno linee
h2 = plot(NaN, NaN, 'bx'); % Dummy plot per alpha_L
h3 = plot(NaN, NaN, 'bx'); % Dummy plot per alpha_Y

% Legenda con i valori di alpha_L e alpha_Y senza linea
legend([h1, h2, h3], {'Eigenvalues', ...
    ['\alpha_L = ', num2str(alpha_L)], ...
    ['\alpha_Y = ', num2str(alpha_Y)]}, 'Location', 'Best');

hold off;


%% Calcolo autovalori H2
eig_CT_H2 = eig(A+B*K_c_CT_H2);
spectral_abscissa_H2 = max(eig_CT_H2);   % Spectral Abscissa

%% Plot Autovalori H2
figure;
hold on;
grid on;
axis equal;

% Impostiamo limiti degli assi in base agli autovalori
y_limit = max(abs(imag(eig_CT_H2))) + 0.2;
x_limit = y_limit;

% Limiti degli assi (estesi di un po' rispetto agli autovalori)
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]);

% Disegnare solo gli autovalori
h1 = plot(real(eig_CT_H2), imag(eig_CT_H2), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Autovalori in blu

% Disegnare gli assi
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1); % Asse X
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1); % Asse Y

% Titolo e etichette degli assi
title('LMI H2');
xlabel('Re');
ylabel('Im');

% Legenda per gli autovalori
legend(h1, {'Eigenvalues'}, 'Location', 'Best');

hold off;

%% Calcolo autovalori Mixed
eig_CT_Mixed = eig(A+B*K_c_CT_Mixed);
spectral_abscissa_Mixed = max(eig_CT_Mixed);   % Spectral Abscissa

%% Plot Autovalori Mixed
% Find the maximum absolute real value for symmetric axes
x_limit = max(abs(real(eig_CT_Mixed))) + 0.2;
y_limit = x_limit; % Make the plot square

% Create the figure
figure;
hold on;
grid on;
axis equal;

% Set symmetric axis limits
xlim([-x_limit, x_limit]);
ylim([-y_limit, y_limit]); 

% Define the sector lines
x_mixed = linspace(-x_limit, 0, 100); % Generate points for the lines
y_sector1 = tan(angle) * x_mixed; % First line (+angle w.r.t. x-axis)
y_sector2 = -tan(angle) * x_mixed; % Second line (-angle w.r.t. x-axis)

% Plot the sector with the two boundary lines
h1 = plot(x_mixed, y_sector1, 'r-', 'LineWidth', 1.5); % Positive real part
h2 = plot(x_mixed, y_sector2, 'r-', 'LineWidth', 1.5); % Negative real part

% Plot the eigenvalues in continuous time
h3 = plot(real(eig_CT_Mixed), imag(eig_CT_Mixed), 'bx', 'MarkerSize', 10, 'LineWidth', 2);

% Draw the axes with symmetric length
plot([-x_limit, x_limit], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-y_limit, y_limit], 'k', 'LineWidth', 1);
h_line = plot(-alpha * ones(1, 100), linspace(-y_limit, y_limit, 100), 'r-', 'LineWidth', 1.5); % Reetta rossa perpendicolare

% Title and legend with the sector angle
title('LMI Sector + Performance + Effort');
legend([h1, h3], {sprintf('Sector + Performance + Effort'), 'Eigenvalues'}, 'Location', 'Best');
xlabel('Re');
ylabel('Im');

hold off;

%% Simulation data
Tfinal=10;
T=0:0.01:Tfinal;
% Random initial condition 
min_x0 = -10;
max_x0 = 10;
random_number = min_x0 + (max_x0-min_x0) .* rand(1,1);  % Generate a random number between [min_x0, max_x0]
x0 = repmat(random_number,36,1);

k = 0;
Q = eye(size(G,2));
x_c_CT_H2(:,1)= x0;
a = 0.1; % Definisce l'ampiezza del rumore
% CT Simulation 
for t=T
    k=k+1;

    % state computation
    x_c_free(:,k)=expm(A*t)*x0; % No control
    x_c_CT(:,k)=expm((A+B*K_c_CT)*t)*x0;
    x_c_CT_perf(:,k)=expm((A+B*K_c_CT_perf)*t)*x0;
    x_c_CT_sector(:,k)=expm((A+B*K_c_CT_sector)*t)*x0;
    x_c_CT_effort(:,k)=expm((A+B*K_c_CT_effort)*t)*x0;
    w = a * (2 * rand(2*n,1) - 1); % rumore uniforme
    x_c_CT_H2(:,k) = expm((A+B*K_c_CT_H2)*t) * x_c_CT_H2(:,1) + w;
    x_c_CT_Mixed(:,k)=expm((A+B*K_c_CT_Mixed)*t)*x0;

    % control variable
    u_c_CT(:,k) = K_c_CT * x_c_CT(:,k);
    u_c_CT_perf(:,k) = K_c_CT_perf * x_c_CT_perf(:,k);
    u_c_CT_sector(:,k) = K_c_CT_sector * x_c_CT_sector(:,k);
    u_c_CT_effort(:,k) = K_c_CT_effort * x_c_CT_effort(:,k);
    u_c_CT_H2(:,k) = K_c_CT_H2 * x_c_CT_H2(:,k);
    u_c_CT_Mixed(:,k) = K_c_CT_Mixed * x_c_CT_Mixed(:,k);
end

%% Continuous Time figure
% Primo grafico: posizione lungo X
figure
plot(T, x_c_free(1,:),T, x_c_CT(1,:), T, x_c_CT_perf(1,:), T, x_c_CT_sector(1,:),T,x_c_CT_effort(1,:),T,x_c_CT_H2(1,:),T,x_c_CT_Mixed(1,:)) % Position of the first Mass along x direction
title('CT controller Position in X')  
grid on
legend('No control','CT Stability', 'CT Performance', 'CT Sector', 'CT Effort', 'CT H2', 'CT Mixed') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (X)') % Etichetta dell'asse y

% Secondo grafico: posizione lungo Y
figure
plot(T, x_c_free(3,:),T, x_c_CT(3,:), T, x_c_CT_perf(3,:), T, x_c_CT_sector(3,:),T, x_c_CT_effort(3,:),T, x_c_CT_H2(3,:),T, x_c_CT_Mixed(3,:)) % Position of the first Mass along y direction
title('CT controllers Position in Y')
grid on
legend('No control','CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2', 'CT Mixed') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (Y)') % Etichetta dell'asse y

% Terzo grafico: variabile di controllo U lungo X
figure
plot(T, u_c_CT(1,:),T, u_c_CT_perf(1,:),T, u_c_CT_sector(1,:),T, u_c_CT_effort(1,:),T, u_c_CT_H2(1,:),T, u_c_CT_Mixed(1,:)) % variabile di controllo lungo x
title('CT control variable in X')
grid on
legend('CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2', 'CT Mixed') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Control action (U(x))') % Etichetta dell'asse x

figure
plot(T, u_c_CT(2,:),T, u_c_CT_perf(2,:),T, u_c_CT_sector(2,:),T, u_c_CT_effort(2,:),T, u_c_CT_H2(2,:),T, u_c_CT_Mixed(2,:)) % variabile di controllo lungo y
title('CT control variable in Y')
grid on
legend('CT Stability', 'CT Performance', 'CT Sector','CT Effort', 'CT H2', 'CT Mixed') % Aggiunge la legenda
xlabel('Time (s)') % Etichetta dell'asse y
ylabel('Control action (U(y))') % Etichetta dell'asse y

%% Grafico delle posizioni di tutte le masse
% Grafico posizioni x
figure
hold on
for i=1:4:33
    variable = 'x';
    numb = num2str(ceil(i/4));
    z = [variable, numb];
    plot(T,x_c_CT(i,:), 'DisplayName', z )   
end
grid on
legend show
title('All masses X positions')
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (X)') % Etichetta dell'asse y

% Grafico posizioni y
figure
hold on
for i=3:4:35
    variable = 'y';
    numb = num2str(ceil(i/4));
    z = [variable, numb];
    plot(T,x_c_CT(i,:), 'DisplayName', z)    
end
grid on
legend show
title('All masses Y positions')
xlabel('Time (s)') % Etichetta dell'asse x
ylabel('Position (Y)') % Etichetta dell'asse y

%% Simulation data only for open loop
Tfinal=150;
T=0:0.015:Tfinal;
% Random initial condition 
min_x0 = -10;
max_x0 = 10;
x0 = min_x0 + (max_x0-min_x0) .* rand(36,1);  % Generate a random number between [min_x0, max_x0]


k = 0;
Q = eye(size(G,2));
x_c_CT_H2(:,1)= x0;

% CT Simulation 
for t=T
    k=k+1;

    % state computation
    x_c_free(:,k)=expm(A*t)*x0; % No control
end

figure
hold on
grid on
plot(T, x_c_free(2,:)) % variabile di controllo lungo x
plot(T, x_c_free(4,:)) % variabile di controllo lungo y
title('xv1 and yv1 in OPEN loop')