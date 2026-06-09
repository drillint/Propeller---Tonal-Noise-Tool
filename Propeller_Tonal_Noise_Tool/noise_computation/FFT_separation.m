function [f_norm, SPL_load, SPL_thick] = FFT_separation(filename, BPF)

    % Caricamento dati
    T = readtable(filename);
    p_ref = 20e-6;

    t   = T.t_obs;
    p_L = T{:,3}; % Loading Pressure
    p_T = T{:,4}; % Thickness Pressure

    % Parametri di campionamento
    fs = 1/mean(diff(t));
    N  = length(t);

    % Rimozione componente media
    p_L = detrend(p_L,'constant');
    p_T = detrend(p_T,'constant');

    % Finestra Hann
    w = hann(N);

    % FFT
    YL = fft(p_L .* w);
    YT = fft(p_T .* w);

    % Correzione per la finestra
    CG = mean(w);   % Coherent Gain

    % Spettro monolato
    YL = YL(1:floor(N/2)+1);
    YT = YT(1:floor(N/2)+1);

    % Ampiezza RMS per bin
    A_L = abs(YL)/(N*CG);
    A_T = abs(YT)/(N*CG);

    A_L(2:end-1) = 2*A_L(2:end-1);
    A_T(2:end-1) = 2*A_T(2:end-1);

    % Conversione Peak -> RMS
    A_L = A_L/sqrt(2);
    A_T = A_T/sqrt(2);

    % SPL
    SPL_load  = 20*log10(A_L/p_ref);
    SPL_thick = 20*log10(A_T/p_ref);

    % Asse frequenze
    f = (0:floor(N/2))*fs/N;

    % Normalizzazione con BPF
    f_norm = f/BPF;

end