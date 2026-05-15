function t_avg = naca0018_mean_thickness(c)
% NACA0018_MEAN_THICKNESS Calcola lo spessore medio geometrico del profilo NACA 0018.
%
%   t_avg = naca0018_mean_thickness(c)
%
%   INPUT:
%       c     - Corda del profilo
%
%   OUTPUT:
%       t_avg - Spessore medio integrato

    % Spessore relativo NACA 0018 (18%)
    t = 0.18;

    % Discretizzazione della corda adimensionale (0 to 1)
    x = linspace(0, 1, 1000);

    % Legge dello spessore NACA 4-digit
    yt = 5 * t * ( ...
        0.2969 * sqrt(x) ...
      - 0.1260 * x ...
      - 0.3516 * x.^2 ...
      + 0.2843 * x.^3 ...
      - 0.1015 * x.^4 );

    % Spessore totale (profilo simmetrico: sopra + sotto)
    thickness = 2 * yt;

    % Calcolo dello spessore medio integrato
    % Nota: Usiamo la trasposta ' per assicurarci che siano vettori colonna per simps
    t_avg = simps(x', thickness') * c;

end