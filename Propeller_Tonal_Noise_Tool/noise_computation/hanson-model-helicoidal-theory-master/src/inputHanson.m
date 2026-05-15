%{
Script for processing user inputs for hanson function.

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
%% Input parameters for Hanson model (frequency domain)
%% Discretization
if strcmp(prop.method,'cosine')
    prop.r_muNodes = cosspace(prop.hub,prop.rt,prop.mu+1)./prop.rt;
    prop.xc_Nodes  = cosspace(-0.5,0.5,prop.nu+1)';
elseif strcmp(prop.method,'linear')
    prop.r_muNodes = linspace(prop.hub,prop.rt,prop.mu+1)./prop.rt;
    prop.xc_Nodes  = linspace(-0.5,0.5,prop.nu+1)';
else
    warning('Wrong Input! Calculations will continue with <strong>linear</strong> method')
    prop.r_muNodes = linspace(prop.hub,prop.rt,prop.mu+1)./prop.rt;
    prop.xc_Nodes  = linspace(-0.5,0.5,prop.nu+1)';
end
    prop.r_mu = (prop.r_muNodes(1:end-1)+prop.r_muNodes(2:end))/2;
    prop.xc   = (prop.xc_Nodes(1:end-1) + prop.xc_Nodes(2:end))/2;
    prop.dz   = diff(prop.r_muNodes);
%% Blade properties
inputData.bD = csvread([inputFolder,file.bD]);
prop.bD  = interpData(inputData.bD, prop.r_mu);

inputData.tb = csvread([inputFolder,file.tb]);
prop.tb = interpData(inputData.tb, prop.r_mu);

inputData.FA = csvread([inputFolder,file.FA]);
prop.FA = interpData(inputData.FA, prop.r_mu);

inputData.MCA = csvread([inputFolder,file.MCA]);
prop.MCA  = interpData(inputData.MCA, prop.r_mu);

%% operating conditions
op.Mx = op.Mfl + op.Min;                        %Effective Mach number [-]
op.Mt = sqrt(op.Mht^2 - (op.Mx)^2);     %Tip mach number [-]
op.omega = op.Mt*Amb.c0/prop.rt;                %Rotation speed [rad/s]
op.Mr = sqrt((op.Mx)^2 + (op.Mt.*prop.r_mu).^2); % blade section Mach number [-]
op.beta = sqrt(1-op.Mfl^2);

%% Observer location

%%%% visual coordinates
if isempty(obs.theta) && isempty(obs.r)
    obs.theta = atan2(obs.R,obs.x);    %observer angle in axial direction [rad]
    obs.r = sqrt(obs.R.^2 + obs.x.^2);    % observer distance from propeller center [-]
elseif isempty(obs.R) && isempty(obs.x)
    obs.R = obs.r.*sin(obs.theta); 
    obs.x = obs.r.*cos(obs.theta);
end

%%%% Retarded system
obs.S0 = sqrt(obs.x.^2 + op.beta^2*obs.R.^2);
obs.xr = 1/op.beta^2*(obs.x + op.Mfl*obs.S0); % x in retarded system
obs.rr = sqrt(obs.xr.^2 + obs.R.^2);  % observer distance from propeller center [-]
obs.thetar = atan2(obs.R,obs.xr);    %observer angle in axial direction [rad]


%% Dependent parameters of the model

dep.df = 1-op.Mfl*cos(obs.thetar); %Doppler factor [-]
dep.omegaD = op.omega./dep.df;

%% Propeller loading
ClCdFlag = ~isempty(file.Cl_span) && ~isempty(file.Cd_span);
TQFlag = ~isempty(file.dTNdr) && ~isempty(file.dQNdr);
if  ClCdFlag && ~TQFlag
    inputData.Cl = [csvread([inputFolder,file.Cl_span]);1,0];
    prop.Cl = interpData(inputData.Cl, prop.r_mu);

    inputData.Cd = [csvread([inputFolder,file.Cd_span]);1,0];
    prop.Cd = interpData(inputData.Cd, prop.r_mu);
    
    out = TQLD(prop.Cl,prop.Cd,Amb.rho,op.Mx*Amb.c0,op.omega,prop.r_mu*prop.rt,prop.bD*2*prop.rt,'LD2TQ');
    prop.dTdz = out(1,:).*prop.rt.*prop.B;
    prop.dQdz = out(2,:).*prop.rt.*prop.B;
    clear out
elseif  TQFlag && ~ClCdFlag
    inputData.dTNdr = [csvread([inputFolder,file.dTNdr]); 1,0];
    prop.dTNdr = interpData(inputData.dTNdr, prop.r_mu);
    prop.dTdz = prop.dTNdr.*prop.rt.*prop.B;
    
    inputData.dQNdr = [csvread([inputFolder,file.dQNdr]); 1,0];
    prop.dQNdr = interpData(inputData.dQNdr, prop.r_mu);
    prop.dQdz = prop.dQNdr.*prop.rt.*prop.B;
    
    out = TQLD(prop.dTNdr,prop.dQNdr,Amb.rho,op.Mx*Amb.c0,op.omega,prop.r_mu*prop.rt,prop.bD*2*prop.rt,'TQ2LD');
    prop.Cl = out(1,:);
    prop.Cd = out(2,:);
    clear out
else
    error('Error. \nPlease provide the file name either for CL and CD (file.Cl_span, file.Cd_span), or for Thrust and Torque (file.dTNdr,file.dQNdr).\n')
end

%% Chordwise Loading
if strcmp(prop.thick_dist,'custom')
    inputData.t       = csvread([inputFolder,file.thick_dist]);
    prop.chord.rR     = inputData.t(1,2:end);
    prop.chord.xc     = inputData.t(2:end,1);
    prop.chord.t_tmax = inputData.t(2:end,2:end);
    prop.HX     = interp2(prop.chord.rR,prop.chord.xc,prop.chord.t_tmax,prop.r_mu,prop.xc,'spline');
end
if strcmp(prop.load_dist,'custom')
    inputData.dCl_dxc   = csvread([inputFolder,file.Cl_chord]);
    prop.chord.dCl.rR   = inputData.dCl_dxc(1,2:end);
    prop.chord.dCl.xc   = inputData.dCl_dxc(2:end,1);
    prop.chord.dCl.data = inputData.dCl_dxc(2:end,2:end);
    prop.fLX            = interp2(prop.chord.dCl.rR, prop.chord.dCl.xc,...
                          prop.chord.dCl.data,prop.r_mu,prop.xc,'spline');
    prop.fLX            = prop.fLX./simps(prop.xc,prop.fLX,1); %Normalize the area to unity
    
    inputData.dCd_dxc   = csvread([inputFolder,file.Cd_chord]);
    prop.chord.dCd.rR   = inputData.dCd_dxc(1,2:end);
    prop.chord.dCd.xc   = inputData.dCd_dxc(2:end,1);
    prop.chord.dCd.data = inputData.dCd_dxc(2:end,2:end);
    prop.fDX            = interp2(prop.chord.dCd.rR, prop.chord.dCd.xc,...
                          prop.chord.dCd.data,prop.r_mu,prop.xc,'spline');
    prop.fDX            = prop.fDX./simps(prop.xc,prop.fDX,1); %Normalize the area to unity
end