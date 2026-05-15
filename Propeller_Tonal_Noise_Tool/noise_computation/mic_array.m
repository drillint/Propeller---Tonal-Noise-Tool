clc; clearvars; close all;

%% Load Data
data{1} = load('4412_BemResults_N900_xtr1.00.mat'); % U=12 m/s, 6000 rpm, 2 blades
data{2} = load('0018_BemResults_N900_xtr1.00.mat'); % U=20m/s, 4000 rpm, 3 blades, sweep -1, chord +0

case_labels = {'NACA 4412', 'NACA 0018'};
colors = {'b', 'r'}; 

%% Global Setup
mic_radius = 1;
mic_number = 100;
angle_step = 2*pi/mic_number;
p_ref = 20e-6;

figure('Color', 'w', 'Position', [100, 100, 1500, 700]);

for dir_mode = 1:2
    ax = subplot(1, 2, dir_mode, polaraxes); 
    hold(ax, 'on');
    
    if dir_mode == 1
        theta_dir = true; % Lateral
    else
        theta_dir = false; % Azimuthal
    end

    all_oaspl_plot = []; 

    for i = 1:2
        bem = data{i}.bem;
        prop = data{i}.prop;
        air = data{i}.air;

        % Parameters
        BPF = bem.omega*prop.B/(2*pi);
        n = 6*prop.B;
        chord = prop.geom.data(4,:) / 1000;
        
        if i == 1
            t_avg = naca4412_mean_thickness(chord);
        else
            t_avg = naca0018_mean_thickness(chord);
        end
        
        c_inf = sqrt(air.gamma*air.R/29*air.Tinf);
        
        % Loop Microfoni
        alfa = 0;
        n_mics = mic_number + 1;
        OASPL_vec = zeros(n_mics, 1);
        angles_vec = zeros(n_mics, 1);

        for m = 1:n_mics
            if theta_dir == true
                theta0 = alfa; phi0 = 0;
            else
                theta0 = pi/2; phi0 = alfa;
            end
            
            phi_c = 0;
            psi = phi_c - phi0 + pi/2;
            
            % Sources
            dL = bem.dT ./ prop.dr;
            dQ = bem.dQ ./ prop.dr;
            dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
            dD = dQ ./ prop.r;

            % ... (inizializzazione parametri)
            k_harmonics = 1:6; % armoniche da risolvere
            harmonics = k_harmonics * prop.B;
            p_rms_bpf = zeros(length(harmonics), 1);

            for idx = 1:length(harmonics)
                j = harmonics(idx); % Solo multipli di B

                % Sorgenti (Thickness + Loading)
                S_prime_j = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.* ...
                    besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf) * exp(1i*j*psi);

                Q_prime_j = dL*cos(theta0).*besselj(j, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*j*psi) ...
                    + dD*sin(theta0)/(2i).*(besselj(j-1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi) ...
                    - besselj(j+1, j*bem.omega.*prop.r*sin(theta0)/c_inf)*exp(1i*(j)*psi));

                I_j = trapz(prop.r, S_prime_j + Q_prime_j);

                %
                C_j = (1i*j*bem.omega*exp(1i*j*bem.omega*mic_radius/c_inf)/(4*pi*mic_radius*c_inf)) * ...
                    I_j * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));

                p_rms_bpf(idx) = abs(C_j)/sqrt(2);
            end

            % Calcolo OASPL delle armoniche calcolate
            SPL_harmonics = 20*log10(p_rms_bpf / p_ref);
            OASPL_val = 10*log10(sum(10.^(SPL_harmonics / 10)));

            OASPL_val(OASPL_val < 0) = 0;
            
            OASPL_vec(m) = OASPL_val;
            angles_vec(m) = alfa;
            alfa = alfa + angle_step;
        end
        
        polarplot(ax, angles_vec, OASPL_vec, 'LineWidth', 2, 'Color', colors{i}, ...
                  'DisplayName', [case_labels{i}, ' OASPL']);
        all_oaspl_plot = [all_oaspl_plot; OASPL_vec];
    end

    min_disp = max(0, floor(min(all_oaspl_plot)/10)*10 - 20);
    max_disp = ceil(max(all_oaspl_plot)/10)*10 + 10;
    rlim(ax, [min_disp, max_disp]);

    r_center = min_disp; 
    
    if theta_dir == true
        ax.ThetaDir = 'counterclockwise';
        ax.ThetaZeroLocation = 'left';   
        polarplot(ax, 0, r_center, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
        title(ax, '\textbf{Lateral Directivity ($\theta$)}', 'Interpreter', 'latex', 'FontSize', 14);
    else
        ax.ThetaDir = 'counterclockwise';
        ax.ThetaZeroLocation = 'top';    
        polarplot(ax, 0, r_center, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
        title(ax, '\textbf{Azimuthal Directivity ($\phi$)}', 'Interpreter', 'latex', 'FontSize', 14);
    end
    grid(ax, 'on');
    if dir_mode == 2, legend(ax, 'Location', 'northeast'); end
end

sgtitle('Comparison NACA 4412 vs 0018: Directivity Patterns', 'FontWeight', 'bold', 'FontSize', 16);