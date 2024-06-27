%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Script to create extra figures
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
format long ;





addpath('../') % add path for function parameters

dim = 2;   % Dimension of the CT (dimension of the BM will be dim-1)
test = 3;  % Test case: healthy = 1, tumor = 2, tumor+SF = 3
Tf = 7;    % Final time (days)
Tp = Tf+1; % Number of time points to store solution
L = 0.1;   % Size of conjunctive tissue (length in dm)
Nx = 40;   % Number of spatial grid points (in each direction if dim == 2)

dx = L/Nx; % The grid cell width
if dim == 1
    x = linspace(0,L,Nx);
elseif dim == 2
    x = linspace(0,L,Nx);
    y = linspace(0,L,Nx);
    [X,Yx] = meshgrid(x,y);
else
    error("Wrong dimension")
end

% Parameter values
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(test,L,dim,Nx,x);

% Initial conditions
IC = initial_conditions(Nx,test,dim);

% Find location in x1 where we find SF
if isempty(find(x==posSF(1)))
    % warning('Given position of SF is NOT a grid point (i.e. cannot have solution at x1=x1SF, we pick the closest x1)')
    x1sf_pos = find(abs(x-posSF(1))==min(abs(x-posSF(1))));
else
    x1sf_pos = find(x==posSF(1));
end
x1sf_pos = x1sf_pos(1);

%% Parameter space to explore for plot
rt_i = linspace(0,rt*100,26);
rp_i = linspace(0,rp*5,26);
x2_SF = [L-posSF(2),L/2.0,posSF(2)];


posSF(2) = x2_SF(1);
BM_far = zeros(length(rt_i),length(rp_i));

for j=1:length(rt_i)
    for k=1:length(rp_i)
        k
        
        % Save
        try
            load(['../Saved_data/Data_extra_plot/far',num2str(j),'_',num2str(k),'.mat'])
        catch
            "data not found... I run simulations"
            % Update secretome source
            [S_t,S_p] = Update_SF_sources_2D(spread_secretum,posSF,rt_i(j),rp_i(k),X,Yx,Nx,sph_t,sph_p);
            % Simulate
            [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
                = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
                rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
                beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
                width_BM );
           
            save(['../Saved_data/Data_extra_plot/far',num2str(j),'_',num2str(k),'.mat'],...
            't_BM','ct_conj','cp_conj','cp_BM','ct_BM','ca_BM','M_BM','cm_BM','cd_BM','c1_BM','c2_BM','c3_BM','ctp_BM','cta_BM');
        
        end 
        % Store for plot
        BM_far(j,k) = M_BM(end,x1sf_pos);
    end
end

posSF(2) = x2_SF(2);
BM_mid = zeros(length(rt_i),length(rp_i));

for j=1:length(rt_i)
    for k=1:length(rp_i)
        % % Update secretome source
        [S_t,S_p] = Update_SF_sources_2D(spread_secretum,posSF,rt_i(j),rp_i(k),X,Yx,Nx,sph_t,sph_p);
        % Simulate
        [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
            = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
            rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
            beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
            width_BM );
        % Save
        save(['../Saved_data/Data_extra_plot/mid',num2str(j),'_',num2str(k),'.mat'],...
            't_BM','ct_conj','cp_conj','cp_BM','ct_BM','ca_BM','M_BM','cm_BM','cd_BM','c1_BM','c2_BM','c3_BM','ctp_BM','cta_BM');
%         load(['../Saved_data/Data_extra_plot/mid',num2str(j),'_',num2str(k),'.mat'])
        % Store for plot
        BM_mid(j,k) = M_BM(end,x1sf_pos);
    end
end

posSF(2) = x2_SF(3);
BM_close = zeros(length(rt_i),length(rp_i));

for j=1:length(rt_i)
    for k=1:length(rp_i)
        % % Update secretome source
        [S_t,S_p] = Update_SF_sources_2D(spread_secretum,posSF,rt_i(j),rp_i(k),X,Yx,Nx,sph_t,sph_p);
        % Simulate
        [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
            = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
            rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
            beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
            width_BM );
        % Save
        save(['../Saved_data/Data_extra_plot/close',num2str(j),'_',num2str(k),'.mat'],...
            't_BM','ct_conj','cp_conj','cp_BM','ct_BM','ca_BM','M_BM','cm_BM','cd_BM','c1_BM','c2_BM','c3_BM','ctp_BM','cta_BM');
        % load(['../Saved_data/Data_extra_plot/close',num2str(j),'_',num2str(k),'.mat'])
        % Store for plot
        BM_close(j,k) = M_BM(end,x1sf_pos);
    end
end

%% Plot

save(['../Saved_data/Data_extra_plot_tf',num2str(Tf),'_Nx',num2str(Nx),'.mat'],'rt_i','rp_i','BM_far','BM_mid','BM_close')

hfig = figure;

% Create a custom colormap where the color corresponding to 0.5*Mmax is black
cmap = parula(256); % Or any other colormap
idx = floor(0.5 * size(cmap, 1)) + 1;
cmap(idx, :) = [0, 0, 0]; % Set color at 0.5*Mmax to black

subplot(1,3,1)
surf(rp_i,rt_i,BM_far)
view(2)
ylim([rt_i(1) rt_i(end)])
xlim([rp_i(1) rp_i(end)])
xlabel('$r_p$ (nM/s)','Interpreter','latex')
ylabel('$r_t$ nM/s)','Interpreter','latex')
caxis([0,Mmax])
colormap(cmap)
colorbar
shading interp
axis square
title('SF at $x_2^{SF}=0.005$ dm')

subplot(1,3,2)
surf(rp_i,rt_i,BM_mid)
view(2)
ylim([rt_i(1) rt_i(end)])
xlim([rp_i(1) rp_i(end)])
xlabel('$r_p$ (nM/s)','Interpreter','latex')
ylabel('$r_t$ nM/s)','Interpreter','latex')
caxis([0,Mmax])
colormap(cmap)
colorbar
shading interp
axis square
title('SF at $x_2^{SF}=0.05$ dm')

subplot(1,3,3)
surf(rp_i,rt_i,BM_close)
view(2)
ylim([rt_i(1) rt_i(end)])
xlim([rp_i(1) rp_i(end)])
xlabel('$r_p$ (nM/s)','Interpreter','latex')
ylabel('$r_t$ nM/s)','Interpreter','latex')
caxis([0,Mmax])
colormap(cmap)
colorbar
ylabel(colorbar, 'BM density at $x_1=x_1^{SF}$ (nM)','Interpreter','latex');
shading interp
axis square
title('SF at $x_2^{SF}=0.095$ dm')


set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','FontSize'),'FontSize',16) 
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')



%% Annexed functions

function [S_t,S_p] = Update_SF_sources_2D(spread_secretum,posSF,rt,rp,X,Yx,Nx,sph_t,sph_p)

    for ii = 1:size(posSF,1)
        %xSF =  reshape(X,[Nx^2,1]);
        xSF = X;
        %ySF =  reshape(Yx,[Nx^2,1]);
        ySF = Yx;
        S_t = (sph_t + rt./(spread_secretum * sqrt(2*pi))*exp(...
            -(xSF-posSF(ii,1)).^2/(2*spread_secretum^2) -(ySF-posSF(ii,2)).^2/(2*spread_secretum^2) )); 
        S_p = (sph_p + rp./(spread_secretum * sqrt(2*pi))*exp(...
        -(xSF-posSF(ii,1)).^2/(2*spread_secretum^2)-(ySF-posSF(ii,2)).^2/(2*spread_secretum^2) )); 
    end
    % recast S_t and S_p in array form
    S_t_arr = zeros(Nx*Nx,1);
    S_p_arr = zeros(Nx*Nx,1);
    for ii = 1:size(S_t,1)
        S_t_arr((ii-1)*Nx+1:(ii)*Nx) = S_t(ii,:);
        S_p_arr((ii-1)*Nx+1:(ii)*Nx) = S_p(ii,:);
    end
    S_t = S_t_arr;
    S_p = S_p_arr;

end

