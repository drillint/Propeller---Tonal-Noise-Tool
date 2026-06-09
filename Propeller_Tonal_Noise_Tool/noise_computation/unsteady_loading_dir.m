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

%% Unsteady Case (k = [-6, 0, 6])
k_vec = [-6, 0, 6];
A = 0.03; 
dL_unsteady = dL0.*A/2; dD_unsteady = dD0.*A/2; 
dL = [dL_unsteady; dL0; dL_unsteady];
dD = [dD_unsteady; dD0; dD_unsteady];

%% Calcolo Direttività Polare OASPL
theta_vec = linspace(0, 2*pi, 361); 
OASPL_steady_theta = zeros(size(theta_vec));
OASPL_unsteady_theta = zeros(size(theta_vec));

for t_idx = 1:length(theta_vec)
    th = theta_vec(t_idx);
    
    S_st = zeros(n, length(prop.r)); Q_st = zeros(n, length(prop.r));
    S_un = zeros(n, length(prop.r)); Q_un = zeros(n, length(k_vec), length(prop.r));
    
    for j = 1:n
        arg_curr = j * bem.omega .* prop.r * sin(th) / c_inf;
        
        S_st(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_curr) * exp(1i*j*psi);
        term_L_st = dL0 .* cos(th) .* besselj(j, arg_curr);
        term_D_st = -dD0 .* sin(th) / (2i) .* (besselj(j-1, arg_curr) - besselj(j+1, arg_curr));
        Q_st(j,:) = (term_L_st + term_D_st) * exp(1i*j*psi);
        
        S_un(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_curr) * exp(1i*j*psi);
        for k_idx=1:1:3
            Q_un(j,k_idx,:) = (dL(k_idx,:)*cos(th).*besselj(j-k_vec(k_idx), arg_curr) ...
                - (dD(k_idx,:)*sin(th)/(2i).*(besselj(j-k_vec(k_idx)-1, arg_curr) ...
                - besselj(j-k_vec(k_idx)+1, arg_curr))))*exp(1i*(j-k_vec(k_idx))*psi);
        end
    end
    
    I_st = trapz(prop.r, (S_st + Q_st), 2);
    I_un = trapz(prop.r, (S_un + squeeze(sum(Q_un, 2))), 2);
    
    C_st_v = zeros(1, n); C_un_v = zeros(1, n);
    for j = 1:n
        blocco_A_curr = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf));
        blocco_D_curr = sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
        C_st_v(j) = blocco_A_curr * I_st(j) * blocco_D_curr;
        C_un_v(j) = blocco_A_curr * I_un(j) * blocco_D_curr;
    end
    
    p_rms_st_v = 2*abs(C_st_v)/sqrt(2);
    p_rms_un_v = 2*abs(C_un_v)/sqrt(2);
    
    val_st = 20*log10(sqrt(sum(p_rms_st_v.^2))/p_ref);
    val_un = 20*log10(sqrt(sum(p_rms_un_v.^2))/p_ref);
    
    if isnan(val_st) || isinf(val_st) || val_st < 0, val_st = 0; end
    if isnan(val_un) || isinf(val_un) || val_un < 0, val_un = 0; end
    
    OASPL_steady_theta(t_idx) = val_st;
    OASPL_unsteady_theta(t_idx) = val_un;
end

%% Calcolo Direttività Azimutale OASPL (360 Gradi)
phi_vec = linspace(0, 2*pi, 361); 
OASPL_steady_phi = zeros(size(phi_vec));
OASPL_unsteady_phi = zeros(size(phi_vec));

theta_plane = pi/2; 
arg_phi_curr = (1:n)' * bem.omega .* prop.r * sin(theta_plane) / c_inf;

for p_idx = 1:length(phi_vec)
    ph = phi_vec(p_idx);
    psi_curr = (pi/2) - ph; 
    
    S_st = zeros(n, length(prop.r)); Q_st = zeros(n, length(prop.r));
    S_un = zeros(n, length(prop.r)); Q_un = zeros(n, length(k_vec), length(prop.r));
    
    for j = 1:n
        S_st(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_phi_curr(j,:)) * exp(1i*j*psi_curr);
        term_L_st = dL0 .* cos(theta_plane) .* besselj(j, arg_phi_curr(j,:));
        term_D_st = -dD0 .* sin(theta_plane) / (2i) .* (besselj(j-1, arg_phi_curr(j,:)) - besselj(j+1, arg_phi_curr(j,:)));
        Q_st(j,:) = (term_L_st + term_D_st) * exp(1i*j*psi_curr);
        
        S_un(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* besselj(j, arg_phi_curr(j,:)) * exp(1i*j*psi_curr);
        for k_idx=1:3
            Q_un(j,k_idx,:) = dL(k_idx,:)*cos(theta_plane).*besselj(j-k_vec(k_idx), arg_phi_curr(j,:))*exp(1i*(j-k_vec(k_idx))*psi_curr) ...
                - (dD(k_idx,:)*sin(theta_plane)/(2i).*(besselj(j-k_vec(k_idx)-1, arg_phi_curr(j,:)) ...
                - besselj(j-k_vec(k_idx)+1, arg_phi_curr(j,:))))*exp(1i*(j-k_vec(k_idx))*psi_curr);
        end
    end
    
    I_st = trapz(prop.r, (S_st + Q_st), 2);
    I_un = trapz(prop.r, (S_un + squeeze(sum(Q_un, 2))), 2);
    
    C_st_v = zeros(1, n); C_un_v = zeros(1, n);
    for j = 1:n
        blocco_A_curr = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf));
        blocco_D_curr = sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
        C_st_v(j) = blocco_A_curr * I_st(j) * blocco_D_curr;
        C_un_v(j) = blocco_A_curr * I_un(j) * blocco_D_curr;
    end
    
    p_rms_st_v = 2*abs(C_st_v)/sqrt(2);
    p_rms_un_v = 2*abs(C_un_v)/sqrt(2);
    
    val_st = 20*log10(sqrt(sum(p_rms_st_v.^2))/p_ref);
    val_un = 20*log10(sqrt(sum(p_rms_un_v.^2))/p_ref);
    
    if isnan(val_st) || isinf(val_st) || val_st < 0, val_st = 0; end
    if isnan(val_un) || isinf(val_un) || val_un < 0, val_un = 0; end
    
    OASPL_steady_phi(p_idx) = val_st;
    OASPL_unsteady_phi(p_idx) = val_un;
end

%% Plot
figure('Color','w','Position', [150, 150, 1300, 550]);
subplot(1,2,1);
polarplot(theta_vec, OASPL_steady_theta, 'b', 'LineWidth', 1.5); hold on;
polarplot(theta_vec, OASPL_unsteady_theta, 'r--', 'LineWidth', 1.5);
ax1 = gca;
ax1.ThetaZeroLocation = 'left';
ax1.ThetaDir = 'counterclockwise';
r_max_t = max([OASPL_steady_theta, OASPL_unsteady_theta]) + 5;
rlim([0, r_max_t]); 
title('\theta Directivity', 'FontSize', 11, 'FontWeight', 'bold');
legend({'Steady', 'Unsteady'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
subplot(1,2,2);
polarplot(phi_vec, OASPL_steady_phi, 'b', 'LineWidth', 1.5); hold on;
polarplot(phi_vec, OASPL_unsteady_phi, 'r--', 'LineWidth', 1.5);
ax2 = gca;
ax2.ThetaZeroLocation = 'top';
ax2.ThetaDir = 'counterclockwise';
r_max_p = max([OASPL_steady_phi, OASPL_unsteady_phi]) + 5;
rlim([0, max(r_max_t, r_max_p)]); 
title('\phi Directivity', 'FontSize', 11, 'FontWeight', 'bold');
legend({'Steady', 'Unsteady'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
sgtitle('OASPL', 'FontSize', 14, 'FontWeight', 'bold');