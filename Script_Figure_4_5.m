%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Thu Oct 23 2023
% Script to generate figures 4 and 5 of the article.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

Nx = 50;
L = 0.1;
dim = 2;
Tf = 7;

%% Simulate test cases and store data for FULL SYSTEM
addpath('./full_system/')  

if dim==1
    addpath('./full_system/1D') 
elseif dim==2
    addpath('./full_system/2D') 
end

% % Simulate healthy test case for full system
try 
    load('Saved_data/Data_test1_2D_Tf7');
catch
    test = 1;
    fun_full_system(dim,test,Nx,Tf,L)
end

try 
    load('Saved_data/Data_test2_2D_Tf7');
catch
    % Simulate tumour test case for full system
    test = 2;
    fun_full_system(dim,test,Nx,Tf,L)
end

try 
    load('Saved_data/Data_test3_2D_Tf7');
catch
    % Simulate tumour+SF test case for full system
    test = 3;
    fun_full_system(dim,test,Nx,Tf,L)
 end


%% PLOTS

% Parameters SF position (case: 1 SF, check consistency with 'parameters.m')
dx = L/Nx; 
x = linspace(dx/2,L-dx/2,Nx);
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(3,L,2,Nx,x);
Mcrit = 0.5*Mmax;

% Find location in x1 where we find SF
if isempty(find(x==posSF(1)))
    % warning('Given position of SF is NOT a grid point (i.e. cannot have solution at x1=x1SF, we pick the closest x1)')
    x1sf_pos = find(abs(x-posSF(1))==min(abs(x-posSF(1))));
else
    x1sf_pos = find(x==posSF(1));
end
x1sf_pos = x1sf_pos(1);

% Load data (2D)
load('Saved_data/Data_test1_2D_Tf7')
M_1 = M_BM(end,:); 
M_1t = M_BM(:,x1sf_pos); 
ct_1 = reshape(ct_conj(end,:),[Nx,Nx]);
cp_1 = reshape(cp_conj(end,:),[Nx,Nx]);
ct_0 = reshape(ct_conj(1,:),[Nx,Nx]);
cp_0 = reshape(cp_conj(1,:),[Nx,Nx]);
load('Saved_data/Data_test2_2D_Tf7');
M_2 = M_BM(end,:); 
M_2t = M_BM(:,x1sf_pos); 
ct_2 = reshape(ct_conj(end,:),[Nx,Nx]);
cp_2 = reshape(cp_conj(end,:),[Nx,Nx]);
load('Saved_data/Data_test3_2D_Tf7');
M_3 = M_BM(end,:); 
M_3t = M_BM(:,x1sf_pos); 
ct_3 = reshape(ct_conj(end,:),[Nx,Nx]);
cp_3 = reshape(cp_conj(end,:),[Nx,Nx]);


%% Plot BM density:  M(t=Tf,x)  and  M(t,x=x_SF)

hfig2 = figure(1);  % save the figure handle in a variable
tunit = 3600*24;

% Along BM length at end time
subplot(1,2,2)
plot(x,M_1,'g','LineWidth',3);
hold on
plot(x,M_2,'b','LineWidth',3);
hold on
plot(x,M_3,'r','LineWidth',3);
hold on
plot([x(1) x(end)], [Mcrit Mcrit], 'k--','LineWidth',2);
text(0, Mcrit, '$M_\mathrm{crit}$', 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','Interpreter','latex');
hold on
plot(posSF(1), 0, 'kx', 'MarkerSize', 10,'LineWidth',3); 
text(posSF(1), 0, '$x_1^{SF}$', 'Interpreter', 'latex', 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
hold on
plot([posSF(1) posSF(1)], [0 Mmax], 'k:','LineWidth',2);
xlabel('$x_1$ (dm)','Interpreter','latex')
ylabel('BM density (nM)')
xlim([0, max(x)])
ylim([0 Mmax])
legend('Healthy','Tumour','Tumour+SF', 'box', 'on', 'Location', 'southwest')
title(['BM at day ',num2str(Tf),' (full system)'])

subplot(1,2,1)
plot(t_BM./(tunit),M_1t,'g','LineWidth',3)
hold on
plot(t_BM./(tunit),M_2t,'b','LineWidth',3)
hold on
plot(t_BM./(tunit),M_3t,'r','LineWidth',3)
hold on
plot([0 t_BM(end)/(tunit)], [Mcrit Mcrit], 'k--','LineWidth',2);
text(0, Mcrit, '$M_\mathrm{crit}$', 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','Interpreter','latex');
xlim([0,t_BM(end)/(tunit)])
ylim([0 Mmax])
xlabel('$t$ (days)','Interpreter','latex')
ylabel('BM density (nM)')
legend('Healthy','Tumour','Tumour+SF', 'box', 'on', 'Location', 'southwest')
title('BM at $x_1=x_1^\mathrm{SF}$ (full system)')

picturewidth = 20; 
hw_ratio = 0.8; 
set(findall(hfig2,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig2,'Units','centimeters','Position',[3 3 2*picturewidth hw_ratio*picturewidth])


%% Plot TIMPs and proMMPs in CT
ctmax = max(max(ct_3));
cpmax = max(max(cp_3));
ctmin = min(min(ct_2));
cpmin = min(min(cp_2));

cpmid = min(min(cp_3));
cpmid = max(max(cp_2));

hfig = figure(2);  

 
% subplot(2,4,1)
% surf(x,x,ct_0'); 
% xlabel('$x_1$')
% ylabel('$x_2$')
% view(0,90)
% shading interp
% axis square
% caxis([ctmin,ctmax])
% title('TIMP-2 at $t=0$h')

% subplot(2,4,2)
subplot(2,3,1)
surf(x,x,ct_1'); 
hold on
plot3(x(:),ones(size(x))*x(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 7.5); % solid line
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
view(2)
shading interp
axis square
axis tight
caxis([ctmin,ctmax])
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['TIMP-2 at $t=$',num2str(Tf*24),'h (healthy)'])
title(['TIMP-2 (healthy)'],'FontSize',24)

% subplot(2,4,3)
subplot(2,3,2)
surf(x,x,ct_2'); 
hold on
pos_dotted = ones(length(x),1).*(rho0==0);
x_solid = x(pos_dotted>0);   
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100, 'r:', 'LineWidth', 7.5); % dotted line
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 7.5); % solid line 
axis tight
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
view(0,90)
shading interp
axis square
caxis([ctmin,ctmax])
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['TIMP-2 at $t=$',num2str(Tf*24),'h (tumour)'])
title(['TIMP-2 (tumour)'],'FontSize',24)

% subplot(2,4,4)
subplot(2,3,3)
surf(x,x,ct_3'); 
hold on
plot3(posSF(1),posSF(2),100, 'k+', 'LineWidth', 3, 'MarkerSize', 8)
% text(posSF(1),posSF(2), 'SF', 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','Interpreter','latex');
hold on
pos_dotted = ones(length(x),1).*(rho0==0);
x_solid = x(pos_dotted>0);   
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100, 'r.', 'LineWidth', 7.5); % dotted line
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 7.5); % solid line 
axis tight
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
view(0,90)
shading interp
axis square
caxis([ctmin,ctmax])
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['TIMP-2 at $t=$',num2str(Tf*24),'h (tumour+SF)'])
title(['TIMP-2 (tumour+SF)'],'FontSize',24)

% subplot(2,4,5)
% surf(x,x,cp_0'); 
% xlabel('$x_1$')
% ylabel('$x_2$')
% view(0,90)
% shading interp
% axis square
% caxis([cpmin,cpmax])
% title('proMMP-2 at $t=0$h')

% subplot(2,4,6)
subplot(2,3,4)
surf(x,x,cp_1'); 
hold on
plot3(x(:),ones(size(x))*x(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 7.5); % solid line
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
view(0,90)
shading interp
axis square
axis tight
caxis([cpmin,cpmid])
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['proMMP-2 at $t=$',num2str(Tf*24),'h (healthy)'])
title(['proMMP-2 (healthy)'],'FontSize',24)

% subplot(2,4,7)
subplot(2,3,5)
surf(x,x,cp_2'); 
hold on
pos_dotted = ones(length(x),1).*(rho0==0);
x_solid = x(pos_dotted>0);   
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100, 'r:', 'LineWidth', 7.5); % dotted line
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 7.5); % solid line 
axis tight
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
view(0,90)
shading interp
axis square
caxis([cpmin,cpmid])
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['proMMP-2 at $t=$',num2str(Tf*24),'h (tumour)'])
title(['proMMP-2 (tumour)'],'FontSize',24)

% subplot(2,4,8)
spend = subplot(2,3,6)
surf(x,x,cp_3'); 
hold on
plot3(posSF(1),posSF(2),100, 'k+', 'LineWidth', 3, 'MarkerSize', 8)
hold on
% plot3(x,ones(size(x))*(x(end)+0.5*dx)-0.0013, ones(size(x))*100, 'r.', 'LineWidth', 6); % dotted line
% plot3(x,ones(size(x))*(x(end)+0.5*dx)-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 6); % solid line 
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100, 'r.', 'LineWidth', 7.5); % dotted line
plot3(x,ones(size(x))*(x(end))-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 7.5); % solid line 
axis tight
view(2)
shading interp
xlabel('$x_1$','FontSize',21)
ylabel('$x_2$','FontSize',21)
axis square
caxis([cpmid,cpmax])
colormap(spend,flipud(autumn))
colorbar
ylabel(colorbar, 'concentration (nM)','Interpreter','latex');
% title(['proMMP-2 at $t=$',num2str(Tf*24),'h (tumour+SF)'])
title(['proMMP-2 (tumour+SF)'],'FontSize',24)


set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','FontSize'),'FontSize',16) 
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
