clc; clearvars; close all;

%% --- LOAD DATA ---
data{1} = load('4412_BemResults_N900_xtr1.00.mat'); 
data{2} = load('0018_BemResults_N900_xtr1.00.mat'); 
case_labels = {'NACA 4412', 'NACA 0018'};
p_ref = 20e-6; 

% Microfono a 90°
theta0 = pi/2; r0 = 1; psi = pi/2;
n_bpf_to_show = 3; % Vogliamo vedere le prime 3 BPF

BPF_ref_4412 = (data{1}.bem.omega/(2*pi)) * data{1}.prop.B;
[f_ref, ref_L, ref_T] = PSD_separation('reference_data_90deg.xlsx', BPF_ref_4412);

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
    dL = bem.dT./prop.dr; dQ = bem.dQ./prop.dr;
    dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
    dD = dQ./prop.r; dD(isnan(dD)) = 0;
    
    % Calcolo solo per j = k*B (armoniche acustiche)
    temp_L = zeros(1, n_bpf_to_show);
    temp_T = zeros(1, n_bpf_to_show);
    
    for k = 1:n_bpf_to_show
        j = k * B; % Solo multipli del numero di pale
        
        arg = j*bem.omega.*prop.r*sin(theta0)/c_inf;
        
        S_prime = 1i*j*bem.omega*air.rho*c_inf*t_avg.*chord.*besselj(j, arg) * exp(1i*j*psi);
        Q_prime = dL*cos(theta0).*besselj(j, arg)*exp(1i*j*psi) + ...
                  dD*sin(theta0)/(2i).*(besselj(j-1, arg)*exp(1i*(j)*psi) - ...
                  besselj(j+1, arg)*exp(1i*(j)*psi));
        
        K_j = (1i*j*bem.omega*exp(1i*j*bem.omega*r0/c_inf)/(4*pi*r0*c_inf)) * sum(exp(-1i*j*(1:prop.B)*2*pi/prop.B));
        
        temp_T(k) = max(0, 20*log10((abs(K_j * trapz(prop.r, S_prime))/sqrt(2))/p_ref));
        temp_L(k) = max(0, 20*log10((abs(K_j * trapz(prop.r, Q_prime))/sqrt(2))/p_ref));
    end
    clean_SPL_L{i} = temp_L;
    clean_SPL_T{i} = temp_T;
end

% --- SUBPLOT 1: VALIDAZIONE (Solo NACA 4412) ---
subplot(2,2,[1 3]); hold on;
plot(f_ref, ref_L, 'r-', 'LineWidth', 1, 'DisplayName', 'Ref: Loading');
plot(f_ref, ref_T, 'b-', 'LineWidth', 1, 'DisplayName', 'Ref: Thickness');

% Plot degli stem
f_bpf = 1:n_bpf_to_show;
stem(f_bpf, clean_SPL_L{1}, 'ro', 'filled', 'LineWidth', 1.2, 'DisplayName', 'An: Loading');
stem(f_bpf, clean_SPL_T{1}, 'bo', 'filled', 'LineWidth', 1.2, 'DisplayName', 'An: Thickness');

title('Validation NACA 4412');
xlabel('f / BPF'); ylabel('SPL [dB]'); grid on;
xlim([0.5, 3.5]); ylim([0, 110]); set(gca, 'XTick', 1:3);
legend('Location', 'northeast');

% --- SUBPLOT 2: CONFRONTO 4412 vs 0018 ---
subplot(2,2,2); hold on;
stem(f_bpf, clean_SPL_L{1}, 'r-o', 'LineWidth', 1.2);
stem(f_bpf, clean_SPL_L{2}, 'r-d', 'LineWidth', 1.2);
title('Loading contribution');
ylabel('SPL [dB]');
xlim([0.5 3.5]); ylim([0 110]);
set(gca,'XTick',1:3); grid on;
legend('4412','0018');

subplot(2,2,4); hold on;
stem(f_bpf, clean_SPL_T{1}, 'b-o', 'LineWidth', 1.2);
stem(f_bpf, clean_SPL_T{2}, 'b-d', 'LineWidth', 1.2);
title('Thickness contribution');
xlabel('f / BPF'); ylabel('SPL [dB]');
xlim([0.5 3.5]); ylim([0 110]);
set(gca,'XTick',1:3); grid on;
legend('4412','0018');

%% Direttività
