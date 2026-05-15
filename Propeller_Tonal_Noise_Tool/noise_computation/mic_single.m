clc; clearvars; close all;

%% Load Data
data{1} = load('4412_BemResults_N900_xtr1.00.mat'); % U=12 m/s, 6000 rpm, 2 blades
data{2} = load('0018_BemResults_N900_xtr1.00.mat'); % U=20m/s, 4000 rpm, 3 blades, sweep -1, chord +0

case_labels = {'NACA 4412', 'NACA 0018'};
colors = {'r', 'b'}; 

p_ref = 20e-6; 
mic = [1, 0, 0; 1, pi/2, 0]; 

BPF_ref = (data{1}.bem.omega/(2*pi)) * data{1}.prop.B;
[f_norm_val, SPL_val_0, SPL_val_90] = PSD_spectrum(BPF_ref);

figure('Color','w','Position', [100, 100, 1300, 600]);


for i = 1:2
    bem = data{i}.bem;
    prop = data{i}.prop;
    air = data{i}.air;

    omega = bem.omega/(2*pi); % [Hz]
    BPF = omega*prop.B;
    n = 6*prop.B; % Numero armoniche rotore
    
    chord = prop.geom.data(4,:); % [mm]
    chord = chord/1000;
    
    if i == 1
        t_avg = naca4412_mean_thickness(chord);
    else
        t_avg = naca0018_mean_thickness(chord);
    end

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
        
        dL = bem.dT./prop.dr;
        dQ = bem.dQ./prop.dr;
        dL(isnan(dL)) = 0;
        dQ(isnan(dQ)) = 0;
        dD = dQ./prop.r;
        
        for j = 1:n
            % Thickness
            S_prime(j,:) = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
                           besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf) * exp(1i*j*psi);
            % Loading
            Q_prime(j,:) = dL*cos(theta0).*besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*j*psi) ...
                + dD*sin(theta0)/(2i).*(besselj(j-1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi) ...
                - besselj(j+1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi));
        end
        
        I_curr = trapz(prop.r, (S_prime + Q_prime), 2)';
        C_curr = zeros(1, n);
        for j = 1:n
            C_curr(j) = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * ...
                        I_curr(j) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
        end
        
        p_rms = abs(C_curr)/sqrt(2);
        SPL_analitico(m, :) = 20*log10(p_rms/p_ref);
    end

    % --- Plot ---
    f_norm_analitica = (1:n) / prop.B; 
    for m = 1:2
        subplot(1, 2, m); hold on;
        if i == 1
            if m == 1, p_ref_plt = plot(f_norm_val, SPL_val_0, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Reference for NACA 4412 (High Fedelity)');
            else, plot(f_norm_val, SPL_val_90, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Reference for NACA 4412 (High Fedelity)'); end
        end
        
        stem(f_norm_analitica, SPL_analitico(m, :), colors{i}, 'filled', 'MarkerSize', 4, 'LineWidth', 1.2, ...
             'DisplayName', [case_labels{i}, ' (Analytic)']);
        
        grid on; grid minor;
        xlabel('f / BPF [-]'); ylabel('SPL [dB]');
        title(['Microfono ', num2str(m)]);
        xlim([0.5, 10]);
        ylim([0, 100]);

        legend('Location', 'northeast');
    end
end
sgtitle('Tonal Noise Code Validation - Comparison NACA 4412 vs 0018', 'FontSize', 14, 'FontWeight', 'bold');