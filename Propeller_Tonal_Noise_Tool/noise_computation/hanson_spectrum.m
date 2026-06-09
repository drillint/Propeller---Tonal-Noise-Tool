 clearvars; clc; close all;

%% Gestione della comunicazione tra il tool Gutin e il tool Hanson 
% 1. Identifichiamo la cartella dove si trova QUESTO file (Hanson_comparison.m)
% mfilename('fullpath') restituisce il percorso del file corrente
progettoPath = fileparts(mfilename('fullpath')); 

% 2. Definiamo il percorso della cartella del modello Hanson
hansonPath = fullfile(progettoPath, 'hanson-model-helicoidal-theory-master');

% 3. Aggiungiamo tutto al path di MATLAB (src, data, ecc.)
addpath(genpath(hansonPath));

% 4. Spostiamoci nella sottocartella 'src' per far trovare i file allo script
cd(fullfile(hansonPath, 'src'));

try
    % 5. Lanciamo il Mainfile di Hanson, fornisce la struttura dati di Output
    % nel workspace, per maggiori informazioni sugli output fare
    % riferimento al readme presente dentro la cartella del tool di Hansojj
    Mainfile;
    
    fprintf('Esecuzione completata con successo.\n');
catch ME
    % Se crasha, torna comunque alla cartella principale prima di dare errore
    cd(progettoPath);
    rethrow(ME);
end

% 6. Torniamo alla cartella principale
% NOTA: Se Mainfile contenesse "clearvars", ricalcoliamo il path al volo
cd(fileparts(mfilename('fullpath')));

% 1. Genera il vettore theta sulla base dell'input in UserInput di Hanson
theta = pi/2;

% 
Loading = SPL.loadm;
Thickness = SPL.thicknessm;
Total = SPL.totalm;

%% CONFRONTO con Tool
% Questa sezione ricalcola per il profilo e lo confronta con la
% struttura SPL calcolata dal tool
data = load('0018_BemResults_N900_xtr1.00.mat'); % U=20m/s, 4000 rpm, 3 blades, sweep -1, chord +0

p_ref = 20e-6;

% Distanza dell'ascoltatore
mic_radius = obs.r*prop.rt; % Coerente con il tuo l'Input in Hanson

% Recupero Input n 
bem_c = data.bem;
prop_c = data.prop;
air_c = data.air;
chord_c = prop_c.geom.data(4,:) / 1000;
t_avg_c = naca0018_mean_thickness(chord_c);
c_inf_c = sqrt(air_c.gamma*air_c.R/air_c.Mw*air_c.Tinf);

% Loop di Calcolo
th0 = theta;
psi = pi/2; % phi0 = 0, phi_c = 0

% Sorgenti
dL = bem_c.dT ./ prop_c.dr / prop.B;
dQ = bem_c.dQ ./ prop_c.dr / prop.B;
dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
dD = dQ ./ prop_c.r;
dL=abs(dL); dD=abs(dD);


k_harmonics = 1:6;
harmonics = k_harmonics * prop_c.B;

p_peak_load  = zeros(length(harmonics), 1);
p_peak_tot  = zeros(length(harmonics), 1);


for idx = 1:length(harmonics)
    j = harmonics(idx);
    arg = j * bem_c.omega .* prop_c.r * sin(th0) / c_inf_c;

    % Sorgenti elementari
    S_j = 1i*j*bem_c.omega*air_c.rho*c_inf_c*t_avg_c.*chord_c.* besselj(j, arg) * exp(1i*j*psi);
    Q_j = dL*cos(th0).*besselj(j, arg)*exp(1i*j*psi) ...
        - dD*sin(th0)/(2i).*(besselj(j-1, arg)*exp(1i*j*psi) - besselj(j+1, arg)*exp(1i*j*psi));

    const = (1i*j*bem_c.omega*exp(1i*j*bem_c.omega*mic_radius/c_inf_c)/(4*pi*mic_radius*c_inf_c)) * sum(exp(-1i*j*(1:prop_c.B)*2*pi/prop_c.B));

    % Calcolo Pressioni di PICCO separate (Convenzione Hanson = 2*abs)
    p_peak_load(idx)  = 2 * abs(const * trapz(prop_c.r, Q_j));
    p_peak_thick(idx) = 2 * abs(const * trapz(prop_c.r, S_j));
    p_peak_tot(idx)   = 2 * abs(const * trapz(prop_c.r, S_j + Q_j));
end

% --- CALCOLO VETTORI SPL 
SPL_load_vec = 20 * log10(p_peak_load / p_ref);
SPL_thick_vec = 20 * log10(p_peak_thick / p_ref);
SPL_tot_vec  = 20 * log10(p_peak_tot / p_ref);

x_axis = 1:length(SPL_tot_vec); 

%% --- PLOT ---
figure('Color', 'w', 'Name', 'SPL Comparison: Tool vs Hanson Formulation', ...
       'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.5]);

comp_names = {'Thickness', 'Loading', 'Thickness + Loading'};
data_hanson = {Thickness(1:length(x_axis)), Loading(1:length(x_axis)), Total(1:length(x_axis))};
data_tool   = {SPL_thick_vec, SPL_load_vec, SPL_tot_vec};
x_shift     = 0.15; 

for s = 1:3
    subplot(1, 3, s); hold on;
    
    stem(x_axis, data_hanson{s}, 'b', 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', 'Hanson Formulation');
    stem(x_axis, data_tool{s}, 'r--', 'Marker', 'x', 'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', 'Tool');
    
    for i = 1:length(x_axis)
        text(x_axis(i) - x_shift, data_hanson{s}(i), sprintf('%.1f', data_hanson{s}(i)), ...
            'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
        text(x_axis(i) + x_shift, data_tool{s}(i), sprintf('%.1f', data_tool{s}(i)), ...
            'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    end
    
    grid on; grid minor; 
    title(comp_names{s});
    xlabel('Harmonic Order (f/BPF) [-]'); 
    ylabel('SPL [dB]');
    xlim([0.5, 3.5]); 
    ylim([0, 70]); 
    set(gca, 'XTick', 1:3); 
    legend('Location', 'northeast');
end