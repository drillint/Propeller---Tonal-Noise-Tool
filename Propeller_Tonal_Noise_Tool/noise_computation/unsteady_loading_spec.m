clc; clearvars; close all;

%% Load Data
data = load('0018_BemResults_N900_xtr1.00.mat'); 

p_ref = 20e-6;

mic = [1, pi/2, 0]; 

bem = data.bem;
prop = data.prop;
air = data.air;

omega_hz = bem.omega/(2*pi);
BPF = omega_hz * prop.B;
n = 6 * prop.B;

chord = prop.geom.data(4,:) / 1000; 
t_avg = naca0018_mean_thickness(chord);
c_inf = sqrt(air.gamma*air.R/29*air.Tinf);

r0 = mic(1); theta0 = mic(2); phi0 = mic(3);
psi = (pi/2) - phi0;

dL0 = bem.dT./prop.dr/prop.B;
dQ0 = bem.dQ./prop.dr/prop.B;
dL0(isnan(dL0)) = 0; dQ0(isnan(dQ0)) = 0;
dD0 = dQ0./prop.r;
dD0(isnan(dD0)) = 0;
dL0=abs(dL0); dD0=abs(dD0);

%% Steady Case (k=0) - theta = pi/2
S_steady = zeros(n, length(prop.r));
Q_steady = zeros(n, length(prop.r));
for j = 1:n
    arg = j * bem.omega .* prop.r * sin(theta0) / c_inf;
    
    S_steady(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
        besselj(j, arg) * exp(1i*j*psi);
    
    term_L = dL0 .* cos(theta0) .* besselj(j, arg);
    term_D = -dD0 .* sin(theta0) / (2i) .* ...
        (besselj(j-1, arg) - besselj(j+1, arg));
        
    Q_steady(j,:) = (term_L + term_D) * exp(1i*j*psi);
end
I_steady = trapz(prop.r, (S_steady + Q_steady), 2);
C_steady = zeros(1, n);
for j = 1:n
    C_steady(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_steady(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
p_rms_steady = 2*abs(C_steady)/sqrt(2);
SPL_steady = 20*log10(p_rms_steady/p_ref);

%% Unsteady Case (k = [-6, 0, 6]) - theta = pi/2
k_vec = [-6, 0, 6];
A = 0.03; 
dL_unsteady = dL0.*A/2; dD_unsteady = dD0.*A/2; 
dL = [dL_unsteady; dL0; dL_unsteady];
dD = [dD_unsteady; dD0; dD_unsteady];

S_prime = zeros(n, length(prop.r));
Q_prime = zeros(n, length(k_vec), length(prop.r));
for j = 1:n
    S_prime(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
        besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf) * exp(1i*j*psi);
    for k_idx=1:3
        Q_prime(j,k_idx,:) = dL(k_idx,:)*cos(theta0).*besselj(j-k_vec(k_idx), j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j-k_vec(k_idx))*psi) ...
            - (dD(k_idx,:)*sin(theta0)/(2i).*(besselj(j-k_vec(k_idx)-1, j*bem.omega.*prop.r*sin(theta0)/c_inf) ...
            - besselj(j-k_vec(k_idx)+1, j*bem.omega.*prop.r*sin(theta0)/c_inf)))*exp(1i*(j-k_vec(k_idx))*psi);
    end
end
Q_total = squeeze(sum(Q_prime, 2));
I_curr = trapz(prop.r, (S_prime + Q_total), 2);
C_curr = zeros(1, n);
for j = 1:n
    C_curr(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_curr(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
p_rms_unsteady = 2*abs(C_curr)/sqrt(2);
SPL_unsteady = 20*log10(p_rms_unsteady/p_ref);

%% Calcolo SPL solo Loading - theta = pi/2
I_steady_load = trapz(prop.r, Q_steady, 2);
C_steady_load = zeros(1, n);
for j = 1:n
    C_steady_load(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_steady_load(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
SPL_steady_load = 20*log10(2*abs(C_steady_load)/sqrt(2)/p_ref);

I_unsteady_load = trapz(prop.r, Q_total, 2);
C_unsteady_load = zeros(1, n);
for j = 1:n
    C_unsteady_load(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
        I_unsteady_load(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
end
SPL_unsteady_load = 20*log10(2*abs(C_unsteady_load)/sqrt(2)/p_ref);

%% Calcolo Spettro ad angolo theta = 0
theta_axis = 0;
S_steady_ax = zeros(n, length(prop.r));
Q_steady_ax = zeros(n, length(prop.r));
S_unsteady_ax = zeros(n, length(prop.r));
Q_unsteady_ax = zeros(n, length(k_vec), length(prop.r));

for j = 1:n
    arg_ax = j * bem.omega .* prop.r * sin(theta_axis) / c_inf;
    
    % Steady axis
    S_steady_ax(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_ax) * exp(1i*j*psi);
    term_L_ax = dL0 .* cos(theta_axis) .* besselj(j, arg_ax);
    term_D_ax = -dD0 .* sin(theta_axis) / (2i) .* (besselj(j-1, arg_ax) - besselj(j+1, arg_ax));
    Q_steady_ax(j,:) = (term_L_ax + term_D_ax) * exp(1i*j*psi);
    
    % Unsteady axis
    S_unsteady_ax(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_ax) * exp(1i*j*psi);
    for k_idx=1:3
        Q_unsteady_ax(j,k_idx,:) = dL(k_idx,:)*cos(theta_axis).*besselj(j-k_vec(k_idx), arg_ax)*exp(1i*(j-k_vec(k_idx))*psi) ...
            - (dD(k_idx,:)*sin(theta_axis)/(2i).*(besselj(j-k_vec(k_idx)-1, arg_ax) ...
            - besselj(j-k_vec(k_idx)+1, arg_ax)))*exp(1i*(j-k_vec(k_idx))*psi);
    end
end

% Integrali e pressioni per caso steady ad angolo theta = 0
I_steady_ax = trapz(prop.r, (S_steady_ax + Q_steady_ax), 2);
I_steady_load_ax = trapz(prop.r, Q_steady_ax, 2);
C_steady_ax = zeros(1, n); C_steady_load_ax = zeros(1, n);

% Integrali e pressioni per caso unsteady ad angolo theta = 0
Q_total_ax = squeeze(sum(Q_unsteady_ax, 2));
I_unsteady_ax = trapz(prop.r, (S_unsteady_ax + Q_total_ax), 2);
I_unsteady_load_ax = trapz(prop.r, Q_total_ax, 2);
C_unsteady_ax = zeros(1, n); C_unsteady_load_ax = zeros(1, n);

for j = 1:n
    blocco_A_ax = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf));
    blocco_D_ax = sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
    
    C_steady_ax(j) = blocco_A_ax * I_steady_ax(j) * blocco_D_ax;
    C_steady_load_ax(j) = blocco_A_ax * I_steady_load_ax(j) * blocco_D_ax;
    
    C_unsteady_ax(j) = blocco_A_ax * I_unsteady_ax(j) * blocco_D_ax;
    C_unsteady_load_ax(j) = blocco_A_ax * I_unsteady_load_ax(j) * blocco_D_ax;
end

SPL_steady_ax = 20*log10(2*abs(C_steady_ax)/sqrt(2)/p_ref);
SPL_steady_load_ax = 20*log10(2*abs(C_steady_load_ax)/sqrt(2)/p_ref);
SPL_unsteady_ax = 20*log10(2*abs(C_unsteady_ax)/sqrt(2)/p_ref);
SPL_unsteady_load_ax = 20*log10(2*abs(C_unsteady_load_ax)/sqrt(2)/p_ref);

f_norm = (1:n) / prop.B;

%% Plot Spettro Acustico (theta = pi/2)
figure('Color','w','Position', [100, 100, 1300, 500]);
plots_data = { {SPL_steady_load, SPL_unsteady_load, 'Loading Noise Only'}, ...
               {SPL_steady, SPL_unsteady, 'Total Noise (Thickness + Loading)'} };

for p = 1:2
    subplot(1,2,p); hold on;
    d = plots_data{p};
    s1 = stem(f_norm, d{1}, 'b', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'b', 'BaseValue', 0, 'DisplayName', 'Steady');
    s2 = stem(f_norm, d{2}, 'r--', 'LineWidth', 1.2, 'Marker', 'x', 'BaseValue', 0, 'DisplayName', 'Unsteady');
    
    for k = 1:2
        idx = d{k} > 1;
        if k == 1
            x_shift = -0.1; align_opt = 'right'; col = 'b';
        else
            x_shift = 0.1; align_opt = 'left'; col = 'r';
        end
        
        text(f_norm(idx) + x_shift, d{k}(idx) + 1.5, ...
             cellstr(num2str(d{k}(idx)', '%.1f')), ...
             'Color', col, 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', align_opt);
    end
    grid on; title(d{3}); legend([s1, s2]); xlim([0.5, 6.5]); 
    ylim([0, max([max(d{1}), max(d{2})])+15]); xlabel('f / BPF [-]'); ylabel('SPL [dB]');
end
sgtitle('Far Field Spectrum (\theta = \pi/2)', 'FontWeight', 'bold');

%% Plot Spettro Acustico (theta = 0)
figure('Color','w','Position', [120, 120, 1300, 500]);
plots_data_ax = { {SPL_steady_load_ax, SPL_unsteady_load_ax, 'Loading Noise Only'}, ...
                  {SPL_steady_ax, SPL_unsteady_ax, 'Total Noise (Thickness + Loading)'} };

for p = 1:2
    subplot(1,2,p); hold on;
    d = plots_data_ax{p};
    s1 = stem(f_norm, d{1}, 'b', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'b', 'BaseValue', 0, 'DisplayName', 'Steady');
    s2 = stem(f_norm, d{2}, 'r--', 'LineWidth', 1.2, 'Marker', 'x', 'BaseValue', 0, 'DisplayName', 'Unsteady');
    
    for k = 1:2
        idx = d{k} > 1;
        if k == 1
            x_shift = -0.1; align_opt = 'right'; col = 'b';
        else
            x_shift = 0.1; align_opt = 'left'; col = 'r';
        end
        
        text(f_norm(idx) + x_shift, d{k}(idx) + 1.5, ...
             cellstr(num2str(d{k}(idx)', '%.1f')), ...
             'Color', col, 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', align_opt);
    end
    grid on; title(d{3}); legend([s1, s2]); xlim([0.5, 6.5]); 
    ylim([0, max([max(d{1}), max(d{2})])+15]); xlabel('f / BPF [-]'); ylabel('SPL [dB]');
end
sgtitle('Far Field Spectrum on Propeller Axis (\theta = 0)', 'FontWeight', 'bold');