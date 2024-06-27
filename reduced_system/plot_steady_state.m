%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Script to Plotting analytic steady state against secretome sources 
% from the CT. 
% It replicates the figures 3 and S1 in the article. 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')


% Physiological levels and sources
ct_ph = 9.6;  
cp_ph = 18.0;
cd_ph = 0.14;

Cd = linspace(0,10*cd_ph,20);
Sp = linspace(0,2*cp_ph,20);
St = linspace(0,1.01*ct_ph,20);

%% Steady state

for i=1:length(St)
    for j=1:length(Sp)
        [ct_ij,cp_ij,ca_ij] = steadystate(Sp(j),St(i),cd_ph);
        ct(i,j) = ct_ij;
        cp(i,j) = cp_ij;
        ca(i,j) = ca_ij;
    end
end
figure(1)
subplot(1,3,1)
surf(Sp,St,ct)
view(0,90)
shading flat
colorbar
axis square
ylabel(colorbar, '$c_t^*$ (nM)','Interpreter','latex');
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('Equilibrium TIMP-2')
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')

subplot(1,3,2)
surf(Sp,St,cp)
view(0,90)
shading flat
colorbar
axis square
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('Equilibrium proMMP-2')
ylabel(colorbar, '$c_p^*$ (nM)','Interpreter','latex');
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')

subplot(1,3,3)
surf(Sp,St,ca)
view(0,90)
shading flat
colorbar
axis square
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('Equilibrium active MMP-2')
ylabel(colorbar, '$c_a^*$ (nM)','Interpreter','latex');
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')
drawnow


%% Steady state: plot varying cd

for i=1:length(St)
    for j=1:length(Sp)
        [ct_ij,cp_ij,ca_ij] = steadystate(Sp(j),St(i),cd_ph*0.5);
        ca1(i,j) = ca_ij;
        [ct_ij,cp_ij,ca_ij] = steadystate(Sp(j),St(i),cd_ph*0.75);
        ca2(i,j) = ca_ij;
    end
end

camax = max(max(ca));


figure(2)
subplot(1,3,1)
surf(Sp,St,ca1)
view(0,90)
shading flat
colorbar
caxis([0, camax]);
axis square
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('$\bar c_d = 50\% \, c_d^{ph}$')
ylabel(colorbar, '$c_a^*$ (nM)','Interpreter','latex');
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')

subplot(1,3,2)
surf(Sp,St,ca2)
view(0,90)
shading flat
colorbar
caxis([0, camax]);
axis square
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('$\bar c_d = 75\% \,  c_d^{ph}$')
ylabel(colorbar, '$c_a^*$ (nM)','Interpreter','latex');
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')

subplot(1,3,3)
surf(Sp,St,ca)
view(0,90)
shading flat
colorbar
caxis([0, camax]);
axis square
ylim([St(1),St(end)])
xlim([Sp(1),Sp(end)])
title('$\bar c_d = 100\%\, c_d^{ph}$')
ylabel(colorbar, '$c_a^*$ (nM)','Interpreter','latex');
ylabel('$\overline{c_t}^*$ (nM)')
xlabel('$\overline{c_p}^*$ (nM)')
drawnow


%% FUNCTIONS

function [ct,cp,ca] = steadystate(sp,st,cd)

% Parameter values
k1 = 2.71;
k2 = 0.14;
k3 = 0.02;
km1 = 1e-4;
km2 = 1e-4;
k123 = k1*k3*cd/(km1*(km2+k3));
KR = k2 + k123*k2;
D_t = 1.29e-6; 
D_p = 1.29e-6; 
D_t_BM = D_t*1e-6; 
D_p_BM = D_p*1e-6;
width_BM = 2e-6; 
kappa_t = 3*D_t_BM/width_BM; 
kappa_p = 3*D_p_BM/width_BM;
hat_kappa_t = kappa_t/width_BM; 
hat_kappa_p = kappa_p/width_BM; 

% Steady state compute
if st==0
    ct = 0;
    cp = sp;
    ca = 0; % though could be anything
else
    amp = (st*hat_kappa_t-sp*hat_kappa_p)/(2*hat_kappa_p);
    kkp = hat_kappa_p/(2*KR);
    kkt = hat_kappa_t/(2*KR);
    ssp = st*(hat_kappa_t)/(KR);
    
    cp = sqrt( ( amp + kkt )^2 + ssp ) - ( amp + kkp );
    ct = st*hat_kappa_t/(hat_kappa_t + KR*cp);
    ca = k123*cp;
end

end

function rup = BM_rupture(sp,st,cd)

    [ct,cp,ca] = steadystate(sp,st,cd);

    % Parameter values
    gamma = 0.236; 
    Km = 0.1357;
    rM = 6.18e-4;
    Mmax = 62.5e3;
    Mcrit = 0.5*Mmax;

    % Condition for BM rupture: rup>0
    rup = ca -  (rM/gamma)*(1+Km/Mcrit)*(1-Mcrit/Mmax);

end

