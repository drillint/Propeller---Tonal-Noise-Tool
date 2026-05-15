%{
Script for calculating frequency domain source function for the thickness.

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
%%
function out = PsiV(kx,varargin)
%% 
% three options for distribution function -
% 1. uniform
% 2. parabolic (default)
% 3. custom

if ~isempty(varargin)
    profile =varargin{1};
    if length(varargin)>1
        X = varargin{2};
        HX = varargin{3};
    end
else
    profile = 'parabolic';
end 

if strcmp(profile,'custom')
   fx = HX.*exp(1i*kx.*X);
   out = simps(X,fx,1);
elseif strcmp(profile,'uniform')
    fx = @(x) exp(1i*kx.*x);
    out = integral(fx,-0.5,0.5,'ArrayValued',true);
else
    fx = @(x) (1-4*x.^2).*exp(1i*kx.*x); 
    out = integral(fx,-0.5,0.5,'ArrayValued',true);
end

