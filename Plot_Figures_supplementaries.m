%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Thu Oct 23 2023
% Script to generate the supplementary figures.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all
close all

Nx = 50;
L = 0.1;
dim = 2;

Figure_to_replicate = "S3" % choices "S3", "S4", "S5" 

if Figure_to_replicate == "S3"
    test = 1;
    try
        load('./Saved_data/Data_test1_2D_Tf7.mat')
    catch
        error("You need to first run the script 'Script_Figure_4_5.m'" + ...
            " with the same parameters as this file")
    end
elseif  Figure_to_replicate == "S4"
    test = 2;
    try
        load('./Saved_data/Data_test2_2D_Tf7.mat')
    catch
        error("You need to first run the script 'Script_Figure_4_5.m'" + ...
            " with the same parameters as this file")
    end
else 
    test = 3;
    try
        load('./Saved_data/Data_test3_2D_Tf7.mat')
    catch
        error("You need to first run the script 'Script_Figure_4_5.m'" + ...
            " with the same parameters as this file")
    end
end

%% Parameters SF position (case: 1 SF, check consistency with 'parameters.m')
dx = L/Nx; 
x = linspace(dx/2,L-dx/2,Nx);
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(3,L,2,Nx,x);
Mcrit = 0.6*Mmax;

% Find location in x1 where we find SF
if isempty(find(x==posSF(1)))
    % warning('Given position of SF is NOT a grid point (i.e. cannot have solution at x1=x1SF, we pick the closest x1)')
    x1sf_pos = find(abs(x-posSF(1))==min(abs(x-posSF(1))));
else
    x1sf_pos = find(x==posSF(1));
end
x1sf_pos = x1sf_pos(1);

tunit = 60*60*24; % Time unit for plotting: 1day

if test == 1
    test_case_str = "healthy";
    test_color = 'g';
elseif test ==2
    test_case_str = "tumor";
    test_color = 'b';
elseif test == 3
    test_case_str = "tumor+SF";
    test_color = 'r';
end

%% Main supplementary plot

t_in = 5;

hfig11 = figure(1);  % save the figure handle in a variable

% TIMP
subplot(4,3,1)
plot(x,ct_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,ct_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,ct_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
% legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
legend('$t=0$', ['$t=$',num2str(round(t_BM(t_in)*24./(tunit),2)),' h'], ['$t=$',num2str(t_BM(100)./tunit),' days'],'SF position','Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['TIMP-2 in the BM'])

% prommp
subplot(4,3,2)
plot(x,cp_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cp_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,cp_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['proMMP-2 in the BM'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% mmp
subplot(4,3,3)
plot(x,ca_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,ca_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,ca_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['Active MMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% mt1-mmp monomer
subplot(4,3,4)
plot(x,cm_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cm_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,cm_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['MT1-MMP monomers'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% mt1-mmp dimer
subplot(4,3,5)
plot(x,cd_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cd_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,cd_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['MT1-MMP dimers'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% complex 1 2MT1-TIMP
subplot(4,3,6)
if test==1
    c1_BM = c1_BM./(1e-9);
end
plot(x,c1_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,c1_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,c1_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
if test==1
    ylabel('Concentration (nM) $\times 10^{-9}$','Interpreter','latex') 
else
    ylabel('Concentration (nM)') 
end
% title(['2MT1-MMP/TIMP-2 complex'])
title(['2MT1-MMP/TIMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% complex 2 2MT1-TIMP-prommp
subplot(4,3,7)
if test==1
    c2_BM = c2_BM./(1e-8);
end
plot(x,c2_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,c2_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,c2_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
if test==1
    ylabel('Concentration (nM) $\times 10^{-8}$','Interpreter','latex') 
else
    ylabel('Concentration (nM)') 
end
% title(['2MT1-MMP/TIMP-2/proMMP-2 complex'])
title(['2MT1-MMP/TIMP-2/proMMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% complex 3 2MT1-TIMP-TIMPS
subplot(4,3,8)
if test==1
    c3_BM = c3_BM./(1e-7);
end
plot(x,c3_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,c3_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,c3_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
if test==1
    ylabel('Concentration (nM) $\times 10^{-7}$','Interpreter','latex') 
else
    ylabel('Concentration (nM)') 
end
% title(['2MT1-MMP/TIMP-2/TIMP-2 complex'])
title(['2MT1-MMP/TIMP-2/TIMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% TIMP-prommp
subplot(4,3,9)
plot(x,ctp_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,ctp_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,ctp_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['TIMP-bound proMMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% TIMP-mmp
subplot(4,3,10)
plot(x,cta_BM(1,:),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cta_BM(t_in,:),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,cta_BM(100,:),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
%legend('$t=0$', ['$t=$',num2str(t_BM(t_in)./tunit),' days'], ['$t=$',num2str(t_BM(100)./tunit),' days'] ,'Interpreter','latex') 
xlabel('$x_1$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['TIMP-bound active MMP-2'])
if test==3
    hold on
    plot(posSF(1), 0, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(1), 0, '$x_1^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% TIMP in CT
subplot(4,3,11)
ct0 = reshape(ct_conj(1,:),[Nx,Nx]);
cti = reshape(ct_conj(t_in,:),[Nx,Nx]);
ctend = reshape(ct_conj(end,:),[Nx,Nx]);
plot(x,ct0(:,x1sf_pos),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cti(:,x1sf_pos),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,ctend(:,x1sf_pos),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
xlabel('$x_2$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['TIMP-2 in CT at $x_1=x_1^{SF}$'])
if test==3
    hold on
    plot(posSF(2)-0.01, 6, 'kx', 'MarkerSize', 6,'LineWidth',2); 
     text(posSF(2)-0.01, 6, '$x_2^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

%proMMP in CT
subplot(4,3,12)
cp0 = reshape(cp_conj(1,:),[Nx,Nx]);
cpi = reshape(cp_conj(t_in,:),[Nx,Nx]);
cpend = reshape(cp_conj(end,:),[Nx,Nx]);
plot(x,cp0(:,x1sf_pos),'Color',test_color,'LineStyle',':','LineWidth',3)
hold on
plot(x,cpi(:,x1sf_pos),'Color',test_color,'LineStyle','--','LineWidth',3)
hold on
plot(x,cpend(:,x1sf_pos),'Color',test_color,'LineStyle','-','LineWidth',3)
axis square
xlabel('$x_2$','Interpreter','latex') 
ylabel('Concentration (nM)') 
title(['proMMP-2 in CT at $x_1=x_1^{SF}$'])
if test==3
    ylim([16 21.5])
    hold on
    plot(posSF(2)-0.01, 16, 'kx', 'MarkerSize', 6,'LineWidth',2); 
    text(posSF(2)-0.01, 16, '$x_2^{SF}$','Interpreter','latex', 'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% ylim([17.6,18.3])




set(findall(hfig11,'-property','FontSize'),'FontSize',10) % adjust fontsize to your document               
set(findall(hfig11,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig11,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
