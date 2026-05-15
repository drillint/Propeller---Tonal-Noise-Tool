close all
clear
clc

%% -- README -- %%
%
% This script generates the polar database for the UI NACA airfoil
% To run you need to put in the same folder (or in the MATLAB path) the function "xfoil.m" and the XFOIL executable "xfoil.exe"
% If you are running on unix, you need to install xfoil from your reference package manager (apt on ubuntu or homebrew on macos, for example)
% Otherwise, you can modify the source code of "xfoil.m" to run the executable (.exe) through wine
%

%% -- Input parameters -- %

% Airfoil of each station
airfoil = 'NACA 0018';

% Angle of attack (AOA) grid in degrees
alpha_min = -5;
alpha_max = 15;
alpha_step = 1;

% Mach number grid
M_min = 0;
M_max = 0.3;
M_step = 0.05;

% Reynolds number grid
Re_min = 1e5;
Re_max = 5e5;
Re_step = 1e5;

%%% -- [DO NOT CHANGE FROM HERE!] -- %%%

% Number of spanwise stations 
N_stations = 25;

% Fixed parameters for xfoil 
N_panels = 220;                         % Number of panels
iterations = 200;                       % Number of iterations
N_crit = 9;                             % N_crit (transition)
xTRU = 1;                               % x of transition on the upper airfoil side (1 = free transition)
xTRL = 1;                               % x of transition on the lower airfoil side (1 = free transition)
xTRflag = (xTRU ~= 1) && (xTRL ~= 1);   % Flag to check if transition is forced on upper OR lower side of the airfoil

% BEMT parameters (they do not enter in xfoil) 
f = 1;                                  % Rotational correction factor
cr2by3 = [1 0.07];                      % Geometric chord parameters
rf123 = [1 1 0];                        % Rotation factors

%% -- Perform calculations [DO NOT CHANGE FROM HERE!] -- %%

% Build AOA grid (we do this for xfoil stability)
alpha_positive = 0:alpha_step:alpha_max;
alpha_negative = 0:-alpha_step:alpha_min;

% Iterate on Re and M 
for Re = Re_min:Re_step:Re_max
    for M = M_min:M_step:M_max
        % Define the output folder name
        output_folder = sprintf('st%d-%d_Re%g-%g_N%g_xTRU%.2f_xTRL%.2f_f%.2f_cr2by3_%.2f-%.2f_rf123_%.2f_%.2f_%.2f_M%.2f',...
                                2, N_stations, Re_min, Re_max, 100*N_crit, xTRU, xTRL, f, cr2by3(1), cr2by3(2), rf123(1), rf123(2), rf123(3), M);
        output_folder = fullfile('bemt','data','PolarData',output_folder);
        
        if Re == Re_min
            % Create the output folder, if not existent; 
            % otherwise, remove the directory and create it again
            if exist(output_folder,'dir')
                rmdir(output_folder, 's')
                mkdir(output_folder)
            else
                mkdir(output_folder)
            end
        end

        % Compute the polar with xfoil (positive AOA values)
        [pol_positive, foil_positive] = xfoil(airfoil, alpha_positive, Re, M, ...
                        sprintf('pane'), ...
                        sprintf('ppar n %d',N_panels), ...
                        sprintf('oper iter %d', iterations), ...
                        sprintf('oper/vpar n %f',N_crit), ...
                        sprintf('oper/vpar xtr %f %f',xTRU, xTRL));

        % Compute the polar with xfoil (negative AOA values)
        [pol_negative, foil_negative] = xfoil(airfoil, alpha_negative, Re, M, ...
                        sprintf('pane'), ...
                        sprintf('ppar n %d',N_panels), ...
                        sprintf('oper iter %d', iterations), ...
                        sprintf('oper/vpar n %f',N_crit), ...
                        sprintf('oper/vpar xtr %f %f',xTRU, xTRL));

        % Combine negative and positive alpha parts to retrieve an ordered
        % alpha array from alpha_min to alpha_max and do the same for the
        % aerodynamic coefficients
        alpha = [flip(pol_negative.alpha); pol_positive.alpha];
        cl = [flip(pol_negative.CL); pol_positive.CL];
        cd = [flip(pol_negative.CD); pol_positive.CD];
        cdp = [flip(pol_negative.CDp); pol_positive.CDp];
        cm = [flip(pol_negative.Cm); pol_positive.Cm];

        % Find and delete duplicates in alpha (mainly alpha = 0)
        [~, unique_alpha_indices] = unique(alpha);
        alpha = alpha(unique_alpha_indices);
        cl = cl(unique_alpha_indices);
        cd = cd(unique_alpha_indices);
        cdp = cdp(unique_alpha_indices);
        cm = cm(unique_alpha_indices);

        for i = 2:N_stations
            % Define the output file name
            output_file_name = sprintf('%s%sst%d_Re%g_N%g_xTRU%.2f_xTRL%.2f_f%.2f_rf123_%.2f_%.2f_%.2f.txt',...
                                        output_folder, filesep, i, Re, 100*N_crit, xTRU, xTRL, f, rf123(1), rf123(2), rf123(3));
    
            % Build header
            header = sprintf('airfoilName %s\nRe %g\nM %g\nN %g\nVACC 0 VACCflag 1\nXTR %g %g XTRflag %d\nalpha,cl,cd,cdp,cm',...
                            airfoil, Re, M, N_crit, xTRU, xTRL, xTRflag);
    
            % Build output matrix
            output_matrix = [alpha, cl, cd, cdp, cm];
    
            % Write the header in the output file
            writematrix(header,output_file_name,'QuoteStrings','none','WriteMode','overwrite')
    
            % Write the output matrix in the output file
            writematrix(output_matrix,output_file_name,'WriteMode','append')
        end
    end
end
