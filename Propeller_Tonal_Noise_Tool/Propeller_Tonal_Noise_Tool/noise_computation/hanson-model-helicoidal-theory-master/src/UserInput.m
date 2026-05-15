%{
Script for setting up the case for Hanson's tonal noise model in frequency domain.

Copyright 2023 Jatinder Goyal

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
%}
%% User Inputs for Hanson model to caculate tonal noise in frequency domain
%%  Ambient conditions
Amb.rho= 1.225;                 % Air density[kg/m^3]
Amb.c0 = 340.2941;              % speed of sound [m/s];
Amb.gamma = 1.4;                % Ratio of specific heats [-]
Amb.pref = 2.e-5;               % Reference pressure for SPL calculations [Pa]

%% Propeller parameters
prop.B = 3;                     % Number of propeller blades [-]
prop.rt = 149.5471e-3;              % Propeller tip radius [m]
prop.hub = 20.886529e-3;               % Propeller hub radius [m]
prop.mu = 25;                  % No. of elements in the radial direction [-]
prop.method = 'linear';         % Two options - 'linear' or 'cosine' for distribution of elements along blade span

%% Blade properties
inputFolder = '../data/Xprop/';   % Folder with blade properties contianing csv files
file.bD = 'bd_rR_0018.csv';          % chord to diamater ratio (csv file)
file.tb = 'tb_rR_0018.csv';          % maximum thickness to chord ratio (csv file)
file.FA = 'FA_rR.csv';          % Face alignment [m] or offset (csv file)
file.MCA = 'MCA_rR.csv';        % midchord alignment or sweep [m] for wept blades (csv file)

prop.thick_dist = 'custom';      % Three options - 'parabolic','uniform','custom' for thickness distrbution along the chord
file.thick_dist = 't_xc_0018.csv'; % Give filename if 'custom' distribution is chosen otherwise you can leave it empty
%% Analysis type - Either provide Cl and Cd files or Thurst and Toqrue files, leaver the other two empty as ''
file.Cl_span = 'Loading/BEM_Cl_rR_0018.csv';   % Lift coefficient distrbution along blade span (without induction)
file.Cd_span = 'Loading/BEM_Cd_rR_0018.csv';   % Drag coefficient distrbution along blade span (without induction)

file.dTNdr = '';        % propeller thrust distribution per blade along radial direction
file.dQNdr = '';        % propeller torque distribution per blade along radial direction

%% chord distribution of forces

prop.load_dist = 'parabolic';      % Three options - 'parabolic','uniform','custom' for loading distrbution along the chord

prop.nu = 1000;                  % No. of elements in the chorddirect   ion [-], relevant for 'custom' only
file.Cl_chord = 'Loading/BEM_dCl_dxc_J1.60_idx71.csv';             % Give filename if 'custom' distribution is chosen otherwise you can leave it empty ''
file.Cd_chord = 'Loading/BEM_dCd_dxc_J1.60_idx71.csv';             % Give filename if 'custom' distribution is chosen otherwise you can leave it empty ''
%% operating condition - Distinction between inflow mach number and flight mach number for Dopller effect
op.Min = 20/Amb.c0;                  % inflow Mach number [-]
op.Mfl = 0;                     % flight Mach number [-] 
op.Mht = sqrt((66.6667*2*pi*prop.rt)^2 + 20^2)/Amb.c0;                 % Helical tip Mach number [-]

%% Observer location - either provide R and X or theta and r, leave other two as empty matrix
obs.R = [];             % prependicualr distance in plane of rotation/ propeller radius [-]
obs.x = [];             % Observer distance in axial direction/ propeller radius [-] 

obs.theta =  (1:1:359).*(pi/180);     % observer angle in axial direction [rad]
obs.r = 1/prop.rt;         % observer distance from propeller center/ propeller radius [-]