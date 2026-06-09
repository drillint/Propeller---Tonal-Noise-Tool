clc; clearvars; close all;

%% --- LOAD DATA ---
data{1} = load('4412_BemResults_N900_xtr1.00.mat'); 
data{2} = load('0018_BemResults_N900_xtr1.00.mat'); 
case_labels = {'NACA 4412', 'NACA 0018'};
p_ref = 20e-6; 

% Microfono a 90°
theta0 = pi/2; r0 = 1; psi = pi/2;
n_bpf_to_show = 4; 

BPF_ref_4412 = (data{1}.bem.omega/(2*pi)) * data{1}.prop.B;
[f_ref, ref_L, ref_T] = FFT_separation('reference_data_90deg.xlsx', BPF_ref_4412);

figure('Color','w','Position', [100, 100, 1400, 650]);

clean_SPL_L = cell(1,2);
clean_SPL_T = cell(1,2);

for i = 1:2
    bem = data{i}.bem; prop = data{i}.prop; air = data{i}.air;
    B = prop.B;
    
    chord = prop.geom.data(4,:)/1000;
    if i == 1, t_avg = naca4412_mean_thickness(chord);
    else, t_avg = naca0018_mean_thickness(chord); end
    
    c_inf = sqrt(air.gamma*air.R/29*air.Tinf);
    dL = bem.dT./prop.dr/B; dQ = bem.dQ./prop.dr/B;
    dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
    dD = dQ./prop.r; dD(isnan(dD)) = 0;
    dL=abs(dL); dD=abs(dD);
    
    % Calcolo solo per j = k*B (armoniche acustiche)
    temp_L = zeros(1, n_bpf_to_show);
    temp_T = zeros(1, n_bpf_to_show);
    
    for k = 1:n_bpf_to_show
        j = k * B; % Solo multipli del numero di pale
        
        arg = j*bem.omega.*prop.r*sin(theta0)/c_inf;
        
        S_prime = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.*besselj(j, arg) * exp(1i*j*psi);
        Q_prime = dL*cos(theta0).*besselj(j, arg)*exp(1i*j*psi) - ...
                  dD*sin(theta0)/(2i).*(besselj(j-1, arg)*exp(1i*(j)*psi) - ...
                  besselj(j+1, arg)*exp(1i*(j)*psi));
        
        K_j = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
        
        temp_T(k) = max(0, 20*log10((2*abs(K_j * trapz(prop.r, S_prime))/sqrt(2))/p_ref));
        temp_L(k) = max(0, 20*log10((2*abs(K_j * trapz(prop.r, Q_prime))/sqrt(2))/p_ref));
    end
    clean_SPL_L{i} = temp_L;
    clean_SPL_T{i} = temp_T;
end
f_bpf = 1:n_bpf_to_show;


% --- SUBPLOT 1: VALIDAZIONE (Solo NACA 4412) ---
figure('Name', 'Validazione NACA 4412', 'Color', 'w');
hold on;
plot(f_ref, ref_L, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Reference: Loading');
plot(f_ref, ref_T, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Reference: Thickness');

f_bpf = 1:n_bpf_to_show;
stem(f_bpf, clean_SPL_L{1}, 'ro', 'filled', 'LineWidth', 1.2, 'DisplayName', 'Tool: Loading');
stem(f_bpf, clean_SPL_T{1}, 'bo', 'filled', 'LineWidth', 1.2, 'DisplayName', 'Tool: Thickness');

idx_L = clean_SPL_L{1} > 0.1; 
text(f_bpf(idx_L) + 0.1, clean_SPL_L{1}(idx_L) + 2.5, ...
    cellstr(num2str(clean_SPL_L{1}(idx_L)', '%.1f')), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'r');

idx_T = clean_SPL_T{1} > 0.1;
text(f_bpf(idx_T) - 0.1, clean_SPL_T{1}(idx_T) - 5.0, ...
    cellstr(num2str(clean_SPL_T{1}(idx_T)', '%.1f')), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');

title('Validation NACA 4412');
xlabel('f / BPF'); ylabel('SPL [dB]'); grid on;
xlim([0.5, 4.5]); ylim([0, 80]); set(gca, 'XTick', 1:3);
legend('Location', 'northeast');


% --- SUBPLOT 2 & 4: CONFRONTO ---
figure('Name', 'Confronto', 'Color', 'w');
for i_plot = 1:2
    subplot(1,2,i_plot); hold on;
    
    if i_plot == 1
        data1 = clean_SPL_L{1}; data2 = clean_SPL_L{2};
        marker1 = 'ro'; marker2 = 'bo'; titolo = 'Loading Noise';
    else
        data1 = clean_SPL_T{1}; data2 = clean_SPL_T{2};
        marker1 = 'rd'; marker2 = 'bd'; titolo = 'Thickness Noise';
    end
    
    stem(f_bpf, data1, marker1, 'LineWidth', 1, 'DisplayName', '4412');
    stem(f_bpf, data2, marker2, 'LineWidth', 1.2, 'DisplayName', '0018');
    
    idx1 = data1 > 0.1; 
    text(f_bpf(idx1) + 0.15, data1(idx1) + 3.0, cellstr(num2str(data1(idx1)', '%.1f')), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'r');
    
    idx2 = data2 > 0.1;
    text(f_bpf(idx2) - 0.15, data2(idx2) - 4.5, cellstr(num2str(data2(idx2)', '%.1f')), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');
    
    title(titolo); xlabel('f / BPF'); ylabel('SPL [dB]');
    xlim([0.5 4.5]); ylim([0 120]); 
    set(gca,'XTick',1:3); grid on;
    legend('4412','0018', 'Location', 'northeast');
end
%% Direttività
clc; clearvars;
data{1} = load('4412_BemResults_N900_xtr1.00.mat'); 
data{2} = load('0018_BemResults_N900_xtr1.00.mat'); 
case_labels = {'NACA 4412', 'NACA 0018'};
colors = {'r', 'b'}; 
styles = {'-', '--'}; % Continuo = Loading, Tratteggiato = Thickness

mic_radius = 1; mic_number = 100; angle_step = 2*pi/mic_number; p_ref = 20e-6;

fig_dir = figure('Color', 'w', 'Name', 'Directivity', 'Position', [100, 100, 1200, 600]);
for dir_mode = 1:2
    theta_dir = (dir_mode == 1);
    ax = subplot(1, 2, dir_mode, polaraxes); hold(ax, 'on');
    
    for i = 1:2
        bem = data{i}.bem; prop = data{i}.prop; air = data{i}.air;
        chord = prop.geom.data(4,:) / 1000;
        t_avg = (i==1)*naca4412_mean_thickness(chord) + (i==2)*naca0018_mean_thickness(chord);
        c_inf = sqrt(air.gamma*air.R/29*air.Tinf);
        
        bem.dT(isnan(bem.dT)) = 0; bem.dQ(isnan(bem.dQ)) = 0;
        dL = bem.dT ./ prop.dr / prop.B; dQ = bem.dQ ./ prop.dr / prop.B; dD = dQ ./ prop.r;
        dL=abs(dL); dD=abs(dD);
        
        OASPL_T = zeros(mic_number + 1, 1); OASPL_L = zeros(mic_number + 1, 1);
        angles_vec = zeros(mic_number + 1, 1);
        
        for m = 1:(mic_number + 1)
            alfa = (m-1) * angle_step;
            if theta_dir, theta0 = alfa; phi0 = 0; else, theta0 = pi/2; phi0 = alfa; end
            psi = - phi0 + pi/2;
            harmonics = (1:6) * prop.B;
            p_rms_T = zeros(length(harmonics), 1); p_rms_L = zeros(length(harmonics), 1);
            
            for idx = 1:length(harmonics)
                j = harmonics(idx); arg = j*bem.omega.*prop.r*sin(theta0)/c_inf;
                S = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.*besselj(j, arg) * exp(1i*j*psi);
                Q = dL*cos(theta0).*besselj(j, arg)*exp(1i*j*psi) - ...
                    dD*sin(theta0)/(2i).*(besselj(j-1, arg)*exp(1i*j*psi) - besselj(j+1, arg)*exp(1i*j*psi));
                K_j = (1i*j*bem.omega*exp(1i*j*bem.omega*mic_radius/c_inf)/(4*pi*mic_radius*c_inf)) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
                p_rms_T(idx) = 2*abs(K_j * trapz(prop.r, S))/sqrt(2);
                p_rms_L(idx) = 2*abs(K_j * trapz(prop.r, Q))/sqrt(2);
            end
            OASPL_T(m) = max(0, 10*log10(sum(10.^((20*log10(p_rms_T/p_ref))/10))));
            OASPL_L(m) = max(0, 10*log10(sum(10.^((20*log10(p_rms_L/p_ref))/10))));
            angles_vec(m) = alfa;
        end
        
        % Plot
        
        polarplot(ax, angles_vec, OASPL_L, [colors{i}, styles{1}], 'LineWidth', 1.2, 'DisplayName', [case_labels{i}, ' Load']);
        polarplot(ax, angles_vec, OASPL_T, [colors{i}, styles{2}], 'LineWidth', 1.2, 'DisplayName', [case_labels{i}, ' Thick']);
    end
    sgtitle('Sources Noise Directivity', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ax.ThetaDir = 'counterclockwise';
    if theta_dir
        ax.ThetaZeroLocation = 'left'; title('Directivity: $\theta$', 'Interpreter', 'latex');
    else
        ax.ThetaZeroLocation = 'top'; title('Directivity: $\phi$', 'Interpreter', 'latex');
    end
    grid(ax, 'on'); legend(ax, 'Location', 'northeastoutside');
end