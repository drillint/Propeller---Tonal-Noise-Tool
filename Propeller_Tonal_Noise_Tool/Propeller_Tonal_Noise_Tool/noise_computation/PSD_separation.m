function [f_norm, SPL_load, SPL_thick] = PSD_separation(filename, BPF)
    % Caricamento dati (con t_obs in Col 1, p_load in Col 3, p_thick in Col 4)
    T = readtable(filename);
    p_ref = 20e-6;
    
    t = T.t_obs;
    p_L = T{:,3}; % Loading Pressure
    p_T = T{:,4}; % Thickness Pressure
    
    % Parametri di campionamento
    fs = 1/mean(diff(t));
    Nw = 1024; % Lunghezza della finestra
    df = fs/Nw;
    
    % Analisi Spettrale (PSD) con metodo di Welch
    % Calcola la PSD, la integraosulla banda df e passala in dB
    [PSD_L, f] = pwelch(detrend(p_L), hann(Nw), Nw/2, Nw, fs);
    SPL_load = 10*log10((PSD_L * df) / p_ref^2);
    
    [PSD_T, ~] = pwelch(detrend(p_T), hann(Nw), Nw/2, Nw, fs);
    SPL_thick = 10*log10((PSD_T * df) / p_ref^2);
    
    % Normalizzazione asse frequenze
    f_norm = f / BPF;
end