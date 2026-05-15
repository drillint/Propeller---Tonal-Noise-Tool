%{
Script for solving Hanson's helicoidal surface theory to estimate the tonal
noise in frequency domain.

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
function P = hanson(Amb,prop,op,obs,dep,m)
nbytes = 0;
for o = 1:length(obs.x)
    for k = 1:length(m)
        var1 = 2*m(k)*prop.B./(op.Mr*dep.df(o));
        var2 = op.Mr.^2.*cos(obs.thetar(o))-op.Mx;
        hm.kx(k,:) = var1.*prop.bD.*op.Mt;
        hm.ky(k,:) = var1.*prop.bD.*var2./prop.r_mu;
        hm.phi0(k,:) = var1.*var2.*prop.FA./(prop.r_mu*2*prop.rt);
        hm.phis(k,:) = var1.*op.Mt.*prop.MCA./(2*prop.rt);
        
        hm.f1(k,:) = -Amb.rho*Amb.c0^2*prop.B*sin(obs.thetar(o)).*exp(1i*m(k)*prop.B.*(dep.omegaD(o)*obs.rr(o)*prop.rt/Amb.c0 - pi/2))...
            ./(8*pi*(obs.R(o)./2)*dep.df(o));
        
        bsl = besselj(m(k)*prop.B,m(k)*prop.B*prop.r_mu*op.Mt*sin(obs.thetar(o))/dep.df(o));
        if strcmp(prop.thick_dist,'custom')
            psiVKx = PsiV(hm.kx(k,:),prop.thick_dist, prop.xc,prop.HX);
        else
            psiVKx = PsiV(hm.kx(k,:),prop.thick_dist);
        end
        if strcmp(prop.load_dist,'custom')
           psiLKx = PsiLoad(hm.kx(k,:),prop.load_dist, prop.xc,prop.fLX); 
           psiDKx = PsiLoad(hm.kx(k,:),prop.load_dist, prop.xc,prop.fDX); 
        else
            psiLoadKx = PsiLoad(hm.kx(k,:),prop.load_dist);
            psiLKx = psiLoadKx;
            psiDKx = psiLoadKx;
        end
        hm.f2(k,:) = op.Mr.^2.*exp(1i*(hm.phi0(k,:) + hm.phis(k,:))).*bsl;
        hm.f3.Vm(k,:) = hm.kx(k,:).^2.*prop.tb.*psiVKx;   
        hm.f3.Dm(k,:) = 1i*hm.kx(k,:).*(prop.Cd/2).*psiDKx;
        hm.f3.Lm(k,:) = -1i*hm.ky(k,:).*(prop.Cl/2).*psiLKx;
        
        prop.phi = asin(op.Mx./op.Mr);
        hm.f4 = prop.B.*prop.bD.*(Amb.rho*Amb.c0^2).*op.Mr.^2.*prop.rt^2;
        hm.f3.Tm(k,:) = -1i./(2*hm.f4).*(hm.kx(k,:).*sin(prop.phi).*psiDKx + hm.ky(k,:).*cos(prop.phi).*psiLKx).*prop.dTdz;
        hm.f3.Qm(k,:) = 1i./(2*hm.f4).*(hm.kx(k,:).*cos(prop.phi).*psiDKx- hm.ky(k,:).*sin(prop.phi).*psiLKx).*prop.dQdz./(prop.rt.*prop.r_mu);

        hm.PVm(k,:) = trapz(prop.r_mu,hm.f1(k,:).*hm.f2(k,:).*hm.f3.Vm(k,:));
        hm.PDm(k,:) = trapz(prop.r_mu,hm.f1(k,:).*hm.f2(k,:).*hm.f3.Dm(k,:));
        hm.PLm(k,:) = trapz(prop.r_mu,hm.f1(k,:).*hm.f2(k,:).*hm.f3.Lm(k,:));
        hm.PTm(k,:) = trapz(prop.r_mu,hm.f1(k,:).*hm.f2(k,:).*hm.f3.Tm(k,:));
        hm.PQm(k,:) = trapz(prop.r_mu,hm.f1(k,:).*hm.f2(k,:).*hm.f3.Qm(k,:));
        %% Print progress
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Progress = %.02f%% \n', 100*(((o-1)*length(m)+k)/(length(m)*length(obs.x))));
        
    end
    P.Vm(:,o) = hm.PVm;
    P.Lm(:,o) = hm.PLm;
    P.Dm(:,o) = hm.PDm;
    P.Tm(:,o) = hm.PTm;
    P.Qm(:,o) = hm.PQm;

end