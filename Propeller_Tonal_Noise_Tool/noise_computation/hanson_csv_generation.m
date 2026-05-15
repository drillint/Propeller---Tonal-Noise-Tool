%% Script per automatizzare la creazione dei CSV per Hanson
% questo file Matlab genera i file csv nel formato richiesto dal solutore
% Hanson di Goyal, specificando databaseNameBEM
% Gli input devono essere coerenti con quelli indicati in UserInput di
% Hanson. I file vengono generati nella cartella hanson_csv_repository, aggiungerli a
% hanson-model-helicoidal-theory-master\data\Xprop per farli leggere al 
% codice sorgente 

databaseNameBEM = '0018_BemResults_N900_xtr1.00.mat'; % deve essere nella stessa cartella, oppure specifica il percorso
%% Prop
DatiBEMT = load(databaseNameBEM); % U=20m/s, 4000 rpm, 3 blades, sweep -1, chord +0
chord = DatiBEMT.prop.geom.data(4,:)/1000;
outputFolder = '../hanson_csv_repository/Prop'; % Cartella di destinazione
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Crea la cartella se non esiste
end

rR = DatiBEMT.prop.rR(:); 


% Formato: {'NomeFile.csv', DatiBEMT.NomeCampo}
mappaFile = {
    'bd_rR_0018.csv', (chord./DatiBEMT.prop.Dia)';
    'tb_rR_0018.csv', 0.18*ones(length(chord),1); % specificare lo spessore del proprio profilo
};

for i = 1:size(mappaFile, 1)
    nomeFile = mappaFile{i, 1};
    valori = mappaFile{i, 2}(:); % Estrae i valori come vettore colonna
    
    % Crea la matrice a due colonne [r/R, Valore]
    matriceCSV = [rR, valori];
    
    % Salva il file CSV
    percorsoCompleto = fullfile(outputFolder, nomeFile);
    writematrix(matriceCSV, percorsoCompleto);
    
    fprintf('File generato: %s\n', percorsoCompleto);
end

% Spessori lungo la corda e lungo il raggio (NACA 0018)

nomeFile = 't_xc_0018.csv'; 

% 1. Parametri
n_c = 25; % stazioni lungo la corda
xc = linspace(-0.5, 0.5, n_c)'; 
x_naca = xc + 0.5; 
T = 0.18; % t/c max

% 2. Calcolo spessore NACA 0018 
t_naca = 5 * T * (0.2969*sqrt(x_naca) - 0.1260*x_naca - 0.3516*x_naca.^2 + ...
                  0.2843*x_naca.^3 - 0.1015*x_naca.^4);

% 3. Costruzione Matrice Dati (t/tmax_globale)
chord_row = chord(:)'; 
rR_row = rR(:)';
n_r = length(rR_row);

% Ogni colonna rappresenta lo spessore fisico t(x) per quella stazione radiale
matrice_fisica = t_naca * chord_row; 

% Calcolo del t_max assoluto di tutta la pala per la normalizzazione
t_max_assoluto = max(matrice_fisica, [], 'all');

matrice_normalizzata = matrice_fisica / t_max_assoluto;
 
% Header: [NaN, rR_1, rR_2, ... rR_n] (1 riga, n_r + 1 colonne)
header = [NaN, rR_row];      

% Corpo: [xc, colonna_1, colonna_2, ... colonna_n] 
corpo = [xc, matrice_normalizzata];      

matrice_finale = [header; corpo];    

% 5. Salvataggio 
percorsoFinale = fullfile(outputFolder, nomeFile); 
writematrix(matrice_finale, percorsoFinale);

%% BEM
outputFolder = '../hanson_csv_repository/BEM'; % Cartella di destinazione

mappaFile = {
    'BEM_Cl_rR_0018.csv', DatiBEMT.bem.Cl';
    'BEM_Cd_rR_0018.csv', DatiBEMT.bem.Cd';
};

for i = 1:size(mappaFile, 1)
    nomeFile = mappaFile{i, 1};
    valori = mappaFile{i, 2}(:);
    
    matriceCSV = [rR, valori];
    
    % Salva il file CSV
    percorsoCompleto = fullfile(outputFolder, nomeFile);
    writematrix(matriceCSV, percorsoCompleto);
    
    fprintf('File generato: %s\n', percorsoCompleto);
end