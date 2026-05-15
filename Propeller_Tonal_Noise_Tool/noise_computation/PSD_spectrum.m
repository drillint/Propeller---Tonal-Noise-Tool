function [f_norm, SPL_0, SPL_90] = PSD_spectrum(BPF)
%% Reference pressure
p_ref = 20e-6; % [Pa]

%% File paths
file0 = 'reference_data_0deg.xlsx';
file90 = 'reference_data_90deg.xlsx';

%% Process Data
files = {file0, file90};

for i = 1:2
    T = readtable(files{i});
    t = T.t_obs;
    p = T.p_tot;
    
    % Sampling
    dt = mean(diff(t));
    fs = 1/dt;
    
    % Pre-processing (Rimozione trend e pulizia dello zero)
    p = detrend(p); 
    
    % Welch parameters
    Nw = 1024;
    df = fs / Nw;
    [PSD, f] = pwelch(p, hann(Nw), Nw/2, Nw, fs);
    
    % SPL definition
    SPL = 10*log10((PSD * df) / p_ref^2);
    SPL_results{i} = SPL; 
    f_norm = f / BPF;
end

% Assegnazione output
SPL_0 = SPL_results{1};
SPL_90 = SPL_results{2};

end