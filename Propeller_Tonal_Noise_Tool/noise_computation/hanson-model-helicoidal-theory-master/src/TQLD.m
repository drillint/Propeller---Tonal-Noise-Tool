%{
Script for processing the loading input from users.

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
function out = TQLD(x1,x2,rho,Uinf,omega,r,b,str)
Vp = sqrt(Uinf^2+(omega*r).^2);
phi = asin(Uinf./Vp);
f = 0.5*rho.*Vp.^2.*b;
if strcmp(str,'TQ2LD')
    dTNdr_f = x1./f;
    dQNdr_fr = x2./(f.*r);
    Cl = cos(phi).*dTNdr_f + sin(phi).*dQNdr_fr;
    Cd = cos(phi).*dQNdr_fr - sin(phi).*dTNdr_f;
    out = [Cl;Cd];
elseif strcmp(str,'LD2TQ')
    Cl = x1; Cd = x2;
    dTNdr_f = cos(phi).*Cl - sin(phi).*Cd;
    dQNdr_fr = sin(phi).*Cl + cos(phi).*Cd;
    out = [dTNdr_f.*f;dQNdr_fr.*(f.*r)];
end