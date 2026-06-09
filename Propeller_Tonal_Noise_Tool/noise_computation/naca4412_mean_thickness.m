function t_avg = naca4412_mean_thickness(c)
% NACA4412_THICKNESS_MEAN
% Calcola lo spessore medio geometrico del profilo NACA 4412
%
% INPUT:
%   c     - corda
%
% OUTPUT:
%   t_avg - spessore medio

% spessore relativo NACA 4412
t = 0.12;

% discretizzazione corda
x = linspace(0,1,1000);

% distribuzione spessore (mezzo profilo NACA 4-digit)
yt = 5*t*( ...
    0.2969*sqrt(x) ...
  - 0.1260*x ...
  - 0.3516*x.^2 ...
  + 0.2843*x.^3 ...
  - 0.1015*x.^4 );

% spessore totale
thickness = 2 * yt;

% integrazione numerica
t_avg = trapz(x, thickness) * c;
end