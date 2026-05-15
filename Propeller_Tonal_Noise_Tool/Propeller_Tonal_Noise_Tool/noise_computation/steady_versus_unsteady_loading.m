clc; clearvars; close all;

% Compare the spectral content at a probe located at 1 meter from the centre in the
% propeller plane with respect to the steady loading configuration at point 2
% Group 3: k=+-6, magnitude of unsteadyness 3%

%% Load Data
data = load('0018_BemResults_N900_xtr1.00.mat'); 
p_ref = 20e-6;
% Probe location
mic = [1, pi/2, 0]; % [r0, theta, phi]

bem = data.bem;
prop = data.prop;
air = data.air;

omega_hz = bem.omega/(2*pi); 
BPF = omega_hz * prop.B;
n = 6 * prop.B; 

chord = prop.geom.data(4,:) / 1000; % [m]
t_avg = naca0018_mean_thickness(chord);
c_inf = sqrt(air.gamma*air.R/29*air.Tinf);

r0 = mic(1); theta0 = mic(2); phi0 = mic(3);
psi = (pi/2) - phi0;

dL0 = bem.dT./prop.dr;
dQ0 = bem.dQ./prop.dr;
dL0(isnan(dL0)) = 0; dQ0(isnan(dQ0)) = 0;
dD0 = dQ0./prop.r;

%% Steady Case (k=0)
% Pre-allocation per il caso stazionario
S_steady = zeros(n, length(prop.r));
Q_steady = zeros(n, length(prop.r));

for j = 1:n
    arg = j * bem.omega .* prop.r * sin(theta0) / c_inf;
    
    % Thickness (Gutin)
    S_steady(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
        besselj(j, arg) * exp(1i*j*psi);
    
    % Loading (Gutin, k=0)
    term_L = dL0 .* cos(theta0) .* besselj(j, arg);
    term_D = dD0 .* sin(theta0) / (2i) .* ...
        (besselj(j-1, arg) - besselj(j+1, arg));
        
    Q_steady(j,:) = (term_L + term_D) * exp(1i*j*psi);
end

I_steady = trapz(prop.r, (S_steady + Q_steady), 2);
C_steady = zeros(1, n);
for j = 1:n
    C_steady(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_steady(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
p_rms_steady = abs(C_steady)/sqrt(2);
SPL_steady = 20*log10(p_rms_steady/p_ref);

%% Unsteady Case (k = [-6, 0, 6])
k_vec = [-6, 0, 6];
A = 0.03; % percentage of unsteadyness w.r.t. radial station
dL_unsteady = dL0.*A/2; dD_unsteady = dD0.*A/2; % value it's halved for the one-siding of the spectrum
dL = [dL_unsteady; dL0; dL_unsteady];
dD = [dD_unsteady; dD0; dD_unsteady];

S_prime = zeros(n, length(prop.r));
Q_prime = zeros(n, length(k_vec), length(prop.r));

for j = 1:n
    % Thickness
    S_prime(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
        besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf) * exp(1i*j*psi);
    % Loading
    for k_idx=1:3
        Q_prime(j,k_idx,:) = dL(k_idx,:)*cos(theta0).*besselj(j-k_vec(k_idx), j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j-k_vec(k_idx))*psi) ...
            + dD(k_idx,:)*sin(theta0)/(2i).*(besselj(j-k_vec(k_idx)-1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j-k_vec(k_idx))*psi) ...
            - besselj(j-k_vec(k_idx)+1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j-k_vec(k_idx))*psi));
    end
end
Q_total = squeeze(sum(Q_prime, 2));
I_curr = trapz(prop.r, (S_prime + Q_total), 2);
C_curr = zeros(1, n);
for j = 1:n
    C_curr(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_curr(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
p_rms_unsteady = abs(C_curr)/sqrt(2);
SPL_unsteady = 20*log10(p_rms_unsteady/p_ref);

%% Calcolo SPL solo Loading
% Caso Steady - Solo Loading
I_steady_load = trapz(prop.r, Q_steady, 2);
C_steady_load = zeros(1, n);
for j = 1:n
    C_steady_load(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_steady_load(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
SPL_steady_load = 20*log10(abs(C_steady_load)/sqrt(2)/p_ref);

% Caso Unsteady - Solo Loading
I_unsteady_load = trapz(prop.r, Q_total, 2);
C_unsteady_load = zeros(1, n);
for j = 1:n
    C_unsteady_load(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_unsteady_load(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
SPL_unsteady_load = 20*log10(abs(C_unsteady_load)/sqrt(2)/p_ref);

%% Plot 
f_norm = (1:n) / prop.B;
figure('Color','w','Position', [100, 100, 1300, 600]);

% Subplot 1: Total Noise (Thickness + Loading)
subplot(1,2,2); hold on;
s1 = stem(f_norm, SPL_steady, 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'BaseValue', 0);
s2 = stem(f_norm, SPL_unsteady, 'r', 'LineWidth', 1.2, 'Marker', 'x', 'BaseValue', 0);
grid on; grid minor;
xlabel('f / BPF [-]'); ylabel('SPL [dB]');
title('Total Noise (Thickness + Loading)');
legend([s1, s2], {'Steady', 'Unsteady'}, 'Location', 'northeast');
xlim([0.5, 6.5]); ylim([0, max(SPL_unsteady) + 10]);

% Subplot 2: Only Loading Noise
subplot(1,2,1); hold on;
s3 = stem(f_norm, SPL_steady_load, 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'BaseValue', 0);
s4 = stem(f_norm, SPL_unsteady_load, 'r', 'LineWidth', 1.2, 'Marker', 'x', 'BaseValue', 0);
grid on; grid minor;
xlabel('f / BPF [-]'); ylabel('SPL [dB]');
title('Loading Noise Only');
legend([s3, s4], {'Steady Loading', 'Unsteady Loading'}, 'Location', 'northeast');
xlim([0.5, 6.5]); ylim([0, max(SPL_unsteady_load) + 10]);

sgtitle('Steady vs Unsteady Comparison (k=\pm6, 3%)', 'FontSize', 14, 'FontWeight', 'bold');