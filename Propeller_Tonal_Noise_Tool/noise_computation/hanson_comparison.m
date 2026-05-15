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
% Questa sezione ricalcola Gutin per il profilo e lo confronta con la
% struttura SPL calcolata dal tool di Hanson
data = load('0018_BemResults_N900_xtr1.00.mat'); % U=20m/s, 4000 rpm, 3 blades, sweep -1, chord +0

p_ref = 20e-6;

% Distanza dell'ascoltatore
mic_radius = obs.r*prop.rt; % Coerente con il tuo l'Input in Hanson

% Recupero Input
bem_c = data.bem;
prop_c = data.prop;
air_c = data.air;
chord_c = prop_c.geom.data(4,:) / 1000;
t_avg_c = naca0018_mean_thickness(chord_c);
c_inf_c = sqrt(air_c.gamma*air_c.R/29*air_c.Tinf);
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
    dL = bem_c.dT ./ prop_c.dr;
    dQ = bem_c.dQ ./ prop_c.dr;
    dL(isnan(dL)) = 0; dQ(isnan(dQ)) = 0;
    dD = dQ ./ prop_c.r;
    
    k_harmonics = 1:6;
    harmonics = k_harmonics * prop_c.B;
    
    p_rms_thick = zeros(length(harmonics), 1);
    p_rms_load  = zeros(length(harmonics), 1);
    p_rms_tot  = zeros(length(harmonics), 1);

    
    for idx = 1:length(harmonics)
        j = harmonics(idx);
        % Argomento Bessel
        arg = j * bem_c.omega .* prop_c.r * sin(th0) / c_inf_c;
        
        % Termine Thickness (S) e Loading (Q)
        S_j = 1i*j*bem_c.omega*air_c.rho*c_inf_c*t_avg_c.*chord_c.* besselj(j, arg) * exp(1i*j*psi);
        Q_j = dL*cos(th0).*besselj(j, arg)*exp(1i*j*psi) ...
            + dD*sin(th0)/(2i).*(besselj(j-1, arg)*exp(1i*j*psi) - besselj(j+1, arg)*exp(1i*j*psi));
        
        % Integrazione
        const = (1i*j*bem_c.omega*exp(1i*j*bem_c.omega*mic_radius/c_inf_c)/(4*pi*mic_radius*c_inf_c)) * sum(exp(-1i*j*(1:prop_c.B)*2*pi/prop_c.B));
        
        p_rms_thick(idx) = abs(const * trapz(prop_c.r, S_j)) / sqrt(2);
        p_rms_load(idx)  = abs(const * trapz(prop_c.r, Q_j)) / sqrt(2);
        p_rms_tot(idx) = abs(const * trapz(prop_c.r, S_j+Q_j)) / sqrt(2);
    end
    
    % Conversione in OASPL [dB]
    OASPL_thickness_calc(m) = 10*log10(sum(10.^((20*log10(p_rms_thick/p_ref))/10)));
    OASPL_loading_calc(m)   = 10*log10(sum(10.^((20*log10(p_rms_load/p_ref))/10)));
    OASPL_total_calc(m)     = 10*log10(sum(10.^((20*log10(p_rms_tot/p_ref))/10)));
end

% Applichiamo il clipping a 0 dB
OASPL_thickness_calc = max(OASPL_thickness_calc, 0);
OASPL_loading_calc   = max(OASPL_loading_calc, 0);
OASPL_total_calc     = max(OASPL_total_calc, 0);

% PLOT DI CONFRONTO (Hanson vs Gutin)
figure('Color', 'w', 'Name', 'Confronto Hanson vs Gutin-Lowson (NACA 0018)');

% --- THICKNESS NOISE ---
subplot(1,3,1, polaraxes);
polarplot(theta, V_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_thickness_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Gutin-Lowson');
title('Thickness Noise'); legend('Location', 'southoutside');

% --- LOADING NOISE ---
subplot(1,3,2, polaraxes);
polarplot(theta, Loading_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_loading_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Gutin-Lowson');
title('Loading Noise'); legend('Location', 'southoutside');

% --- TOTAL NOISE ---
subplot(1,3,3, polaraxes);
polarplot(theta, Total_plot, 'b', 'LineWidth', 1.5, 'DisplayName', 'Hanson'); hold on;
polarplot(theta, OASPL_total_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Gutin-Lowson');
title('Total Noise'); legend('Location', 'southoutside');

sgtitle('Comparison: Hanson Model vs Gutin-Lowson (NACA 0018)', 'FontSize', 16);

%% Volume della pala che vede il metodo Hanson, per distribuzione 'custom'
% Nota: per altre distribuzioni lungo la corda non viene definita prop.HX
% Calcolo del Coefficiente d'Area (ka) per ogni sezione radiale
ka = zeros(1, prop.mu);
for i = 1:prop.mu
    ka(i) = simps(prop.xc, prop.HX(:,i)); % sono le aree adimensionalizzate rispetto alla corda e allo spessore massimo globale
end

% Calcolo delle grandezze fisiche
R_fisico = prop.r_mu * prop.rt;
Corda_fisica = prop.bD * (2 * prop.rt); % prop.bD contiene la corda adimensionalizzata con il diametro

% HX è normalizzata con tmax globale, il valore 1.0 nella matrice rappresenta lo spessore massimo 
% assoluto della pala. Visto che i profili sono 0018 dall' hub al tip,
% allora:
Tmax_assoluto = prop.tb * max(Corda_fisica); % spessore massimo fisico

% L'area fisica si ottiene moltiplicando la corda locale per lo spessore 
% di riferimento, che viene poi scalato dai valori (già ridotti) di ka.
Area_sez = Corda_fisica .* Tmax_assoluto .* ka;

% Volume Totale della singola pala [m^3]
Volume_Pala = simps(R_fisico, Area_sez);

fprintf('Il volume della singola pala per Gutin è: %e m^3\n', Vol_gutin)
fprintf('Il volume della singola pala per Hanson è: %e m^3\n', Volume_Pala)

