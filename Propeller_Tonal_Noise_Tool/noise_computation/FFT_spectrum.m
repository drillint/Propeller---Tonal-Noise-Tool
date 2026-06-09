function [f_norm, SPL_0, SPL_90] = FFT_spectrum(BPF)

%% Reference pressure
p_ref = 20e-6; % [Pa]

%% File paths
file0  = 'reference_data_0deg.xlsx';
file90 = 'reference_data_90deg.xlsx';

%% Process Data
files = {file0, file90};

for i = 1:2

    T = readtable(files{i});

    t = T.t_obs;
    p = T.p_tot;

    %% Sampling
    dt = mean(diff(t));
    fs = 1/dt;

    %% Pre-processing
    p = detrend(p);

    %% FFT parameters
    N = length(p);

    % Hann window
    w = hann(N);

    % Coherent gain correction
    CG = mean(w);

    %% FFT
    Y = fft(p .* w);

    % Single-sided spectrum
    Y = Y(1:floor(N/2)+1);

    % Amplitude spectrum
    A = abs(Y)/(N*CG);
    A(2:end-1) = 2*A(2:end-1);

    % RMS amplitude
    A = A/sqrt(2);

    %% SPL
    SPL = 20*log10(A/p_ref);

    %% Frequency axis
    f = (0:floor(N/2))*fs/N;

    SPL_results{i} = SPL;
    f_norm = f/BPF;

end

%% Outputs
SPL_0  = SPL_results{1};
SPL_90 = SPL_results{2};

end