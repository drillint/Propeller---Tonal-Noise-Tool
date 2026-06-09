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

%% Direttività
% 1. Genera il vettore theta sulla base dell'input in UserInput di Hanson
theta = linspace(0, 2*pi, length(SPL.total));

% 2. Applica il clipping (valori < 0 diventano 0)
V_plot = max(SPL.V, 0);
Loading_plot = max(SPL.loadingLD, 0);
Total_plot = max(SPL.total, 0);

% 3. Plot 
figure('Name', 'OASPL Direttivity');
polarplot(theta, V_plot, 'LineWidth', 1.5); hold on;
polarplot(theta, Loading_plot, 'LineWidth', 1.5);
polarplot(theta, Total_plot, 'k--', 'LineWidth', 2);

title('OASPL Direttivity [dB]');
legend('Thickness', 'Loading', 'Total', 'Location', 'southoutside');

%% CONFRONTO Gutin vs Hanson
% Questa sezione ricalcola iltool per il profilo e lo confronta con la
% struttura SPL calcolata dal tool di Hanson
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
% Volume della pala visto da Gutin
Vol_gutin=simps(prop_c.r,t_avg_c.*chord_c);

% Pre-allocazione vettori
OASPL_thickness_calc = zeros(size(theta));
OASPL_loading_calc   = zeros(size(theta));
OASPL_total_calc     = zeros(size(theta));

% Loop di Calcolo
for m = 1:length(theta)
    th0 = theta(m);
    psi = pi/2; % phi0 = 0, phi_c = 0
    
    % Sorgenti
    dL = bem_c.dT ./ prop_c.dr / prop.B;
    dQ = bem_c.dQ ./ prop_c.dr / prop.B;
    dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
    dD = dQ ./ prop_c.r;
    dL=abs(dL); dD=abs(dD);
    
    k_harmonics = 1:6;
    harmonics = k_harmonics * prop_c.B;
    
    p_peak_thick = zeros(length(harmonics), 1);
    p_peak_load  = zeros(length(harmonics), 1);
    p_peak_tot  = zeros(length(harmonics), 1);

    
    for idx = 1:length(harmonics)
        j = harmonics(idx);
        arg = j * bem_c.omega .* prop_c.r * sin(th0) / c_inf_c;
        
        % Sorgenti elementari
        S_j = 1i*j*bem_c.omega*air_c.rho*c_inf_c*t_avg_c.*chord_c.* besselj(j, arg) * exp(1i*j*psi);
        Q_j = dL*cos(th0).*besselj(j, arg)*exp(1i*j*psi) ...
            - dD*sin(th0)/(2i).*(besselj(j-1, arg)*exp(1i*j*psi) - besselj(j+1, arg)*exp(1i*j*psi));
        
        % Costante comune (Verifica se la formulazione Hanson che stai usando la applica a entrambi!)
        const = (1i*j*bem_c.omega*exp(1i*j*bem_c.omega*mic_radius/c_inf_c)/(4*pi*mic_radius*c_inf_c)) * sum(exp(-1i*j*(1:prop_c.B)*2*pi/prop_c.B));
        
        % Calcolo Pressioni di PICCO separate (Convenzione Hanson = 2*abs)
        p_peak_thick(idx) = 2 * abs(const * trapz(prop_c.r, S_j));
        p_peak_load(idx)  = 2 * abs(const * trapz(prop_c.r, Q_j));
        p_peak_tot(idx)   = 2 * abs(const * trapz(prop_c.r, S_j + Q_j));
    end

    OASPL_thickness_calc(m) = 10*log10(sum(10.^((20*log10(p_peak_thick/p_ref))/10))); 
    OASPL_loading_calc(m)   = 10*log10(sum(10.^((20*log10(p_peak_load/p_ref))/10)));
    OASPL_total_calc(m)     = 10*log10(sum(10.^((20*log10(p_peak_tot/p_ref))/10)));
end

% Applichiamo il clipping a 0 dB
OASPL_thickness_calc = max(OASPL_thickness_calc, 0);
OASPL_loading_calc   = max(OASPL_loading_calc, 0);
OASPL_total_calc     = max(OASPL_total_calc, 0);

% PLOT DI CONFRONTO (Hanson vs tool)
figure('Color', 'w', 'Name', 'Confronto Tool - Hanson (NACA 0018)');

% --- THICKNESS NOISE ---
subplot(1,3,1, polaraxes);
polarplot(theta, V_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_thickness_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Tool');
title('Thickness Noise'); legend('Location', 'southoutside');

% --- LOADING NOISE ---
subplot(1,3,2, polaraxes);
polarplot(theta, Loading_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_loading_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Tool');
title('Loading Noise'); legend('Location', 'southoutside');

% --- TOTAL NOISE ---
subplot(1,3,3, polaraxes);
polarplot(theta, Total_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_total_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Tool');
title('Total Noise'); legend('Location', 'southoutside');

sgtitle('{\theta Directivity}: Tool vs Hanson Formulation - NACA 0018', 'FontSize', 16);

% Applichiamo le modifiche a tutti i subplot
h_axes = findall(gcf, 'type', 'polaraxes');

for i = 1:length(h_axes)
    % Imposta la direzione antioraria (default è 'counterclockwise', ma specifichiamolo)
    set(h_axes(i), 'ThetaDir', 'counterclockwise');
    
    % Imposta lo zero a sinistra (che corrisponde a 180 gradi rispetto allo standard)
    set(h_axes(i), 'ThetaZeroLocation', 'left');
end