clc; clearvars; close all;

%% Load Data
data = load('4412_BemResults_N900_xtr1.00.mat'); 
p_ref = 20e-6; 
mic = [1, 0, 0; 1, pi/2, 0]; 
BPF_ref = (data.bem.omega/(2*pi)) * data.prop.B;
[f_norm_val, SPL_val_0, SPL_val_90] = FFT_spectrum(BPF_ref);

figure('Color','w','Position', [100, 100, 1300, 600]);
bem = data.bem;
prop = data.prop;
air = data.air;
omega = bem.omega/(2*pi); % [Hz]
n = 6*prop.B; 
chord = prop.geom.data(4,:)/1000;
t_avg = naca4412_mean_thickness(chord);
c_inf = sqrt(air.gamma*air.R/29*air.Tinf);

SPL_analitico = zeros(2, n);
for m = 1:2
    r0 = mic(m,1);
    theta0 = mic(m,2);
    phi0 = mic(m,3);
    phi_c = 0;
    psi = phi_c - phi0 + (pi/2);
    
    S_prime = zeros(n, length(prop.r));
    Q_prime = zeros(n, length(prop.r));
    
    dL = bem.dT./prop.dr/prop.B;
    dQ = bem.dQ./prop.dr/prop.B;
    dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0; dD = dQ./prop.r;
    
    for j = 1:n
        S_prime(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
                       besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf) * exp(1i*j*psi);
        Q_prime(j,:) = dL*cos(theta0).*besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*j*psi) ...
            - dD*sin(theta0)/(2i).*(besselj(j-1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi) ...
            - besselj(j+1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi));
    end
    
    I_curr = trapz(prop.r, (S_prime + Q_prime), 2)';
    C_curr = zeros(1, n);
    for j = 1:n
        C_curr(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
                    I_curr(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
    end
    p_rms = 2*abs(C_curr)/sqrt(2);
    SPL_analitico(m, :) = 20*log10(p_rms/p_ref);
end

% --- Plot ---
f_norm_analitica = (1:n) / prop.B; 
for m = 1:2
    subplot(1, 2, m); hold on;
    stem(f_norm_analitica, SPL_analitico(m, :), 'r', 'filled', 'MarkerSize', 4, 'LineWidth', 1.2, ...
         'DisplayName', 'Tool');
    if m == 1
        plot(f_norm_val, SPL_val_0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Mid-high fidelity (Reference)');
    else
        plot(f_norm_val, SPL_val_90, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Mid-high fidelity (Reference)');
    end
    
    if m == 2
        valid_tool = SPL_analitico(m, :) > 10;
        if any(valid_tool)
            pos_x = f_norm_analitica(valid_tool)';
            valori = SPL_analitico(m, valid_tool)';
            text(pos_x + 0.2, valori + 2, cellstr(num2str(valori, '%.1f')), ...
                 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'left');
        end
    end
    
    grid on; grid minor;
    xlabel('f / BPF [-]'); ylabel('SPL [dB]');
    title(['Microphone ', num2str(m)]);
    xlim([0.5, 10]); ylim([0, 100]);
    legend('Location', 'northeast');
end
sgtitle('Tonal Noise Code Validation - NACA 4412', 'FontSize', 14, 'FontWeight', 'bold');sgtitle('Tonal Noise Code Validation - NACA 4412', 'FontSize', 14, 'FontWeight', 'bold');