%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (*) Mathematical modelling of the role of senescent fibroblastsin tumor %
%      initiation and evolution: rupture of the basementmembrane.         %
%                                                                         %
%    Luís Almeida, Alexandre Poulain, Albin Pourtier, Chiara Villa        %
%                                                                         %
% Code copyrights: Alexandre Poulain, Chiara Villa                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Code to conduct global sensitivity analysis on the full and reduced PDE %
% models proposed in the paper above (*), using the SAFE Toolbox          %
% (https://safetoolbox.github.io/) developped by Pianosi et al (**) and   %
% distributed under the GNU Public License Version 3.                     %
%                                                                         %
% (**) Pianosi, F., Sarrazin, F., Wagener, T. (2015), A Matlab toolbox    %
% for Global Sensitivity Analysis, Environmental Modelling & Software,    %
% 70, 80-85. doi.org/10.1016/j.envsoft.2015.04.009                        %
%                                                                         %
% To run this code, download SAFE (Matlab version) from:                  %
% >>> https://github.com/SAFEtoolbox/SAFE-matlab                          %
% and include it in the parent directory of this file                     %
%                                                                         %
% This file also calls the function 'PDE_solve.m' which calls             %
% functions in the folder 'PDE_solver', containing the code to            %
% simulate the full and reduced PDE+ODE systems presented in (*):         %
% please check these files are in the same directory as this file         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% We conduct SA 1) Morris Method (Elementary Effects Test) for initial    %
% screening of most influential parameters.                               %
%                                                                         %
% (***) Pianosi, F., et al. (2016), Sensitivity analysis of               %
% environmental models: A systematic review with practical workflow,      %
% Environmental Modelling & Software 79, 214-232.                         %
% doi.org/10.1016/j.envsoft.2016.02.008                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
format long

%% %%%%%%%%%%%  GENERAL SA SET UP: MODEL AND INPUT  %%%%%%%%%%%%%%%%%%%%%%%
myfun  = 'PDE_solve';


%%% Define output function 
addpath('../');

%%% Add path to the SAFE package
addpath('../../../../SAFE-matlab/sampling');
addpath('../../../../SAFE-matlab/EET'); 
addpath('../../../../SAFE-matlab/visualization');

addpath('../full_system/1D/');

%%% Choose test scenario (healthy, tumour, tumour_SF)
test_case = 3; % heatlthy = 1, tumor = 2, tumor+SF = 3

%%% Choose which PDE system to solve ('full' = 1 or 'reduced' = 2)
system = 1;

%%% Choose parameters to test 
par_case = 2; % 1= reaction rates, 2 = rest of the parameters, 

%%% compute model or load results
compute_mod = 1; % 1 to compute the model


%%% Choose the dimension and fix the spatial discretization
dim = 1;
L = 0.1; % dm
Nx = 20; % Number of spatial grid points (in each direction if dim == 2)
dx = L/Nx; % the grid size
if dim == 1
    x = linspace(0,L,Nx);
elseif dim == 2
    x = linspace(0,L,Nx);
    y = linspace(0,L,Nx);
    [X,Yx] = meshgrid(x,y);
else
    error("Wrong dimension")
end

%%% CHOOSE which set of parameters we are testing

% Baseline parameters
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,pos,rt,rp] = parameters(test_case,L,dim,Nx,x);
%IC = initial_conditions(Nx,test_case,dim);
if system == 2
    cd = 1.e-5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% INPUT A) biochemical reaction parameters, for Morris Method %%%%%%%
if par_case == 1 % the recation rates
    if system == 2 % if reduced system c_d iis an input
        y = [k1, km1, k2, km2, k3, beta_t, beta_p, cd]; % Baseline parameter values
        Y_labels = {'$k_1$','$k_{-1}$',...
           '$k_2$','$k_{-2}$','$k_3$','$\beta_t$','$\beta_p$', '$c_d$' }; % Parameter names
    elseif system == 1
        y = [k0, km0, k1, km1, k2, km2, k3,beta_t,beta_ta,beta_m,beta_d,beta_p,beta_a,beta_tp]; % Baseline parameter values
        Y_labels = {'$k_0$','$k_{-0}$','$k_1$','$k_{-1}$',...
           '$k_2$','$k_{-2}$','$k_3$', '$\beta_t$','$\beta_{ta}$','$\beta_m$','$\beta_d$',...
            '$\beta_p$','$\beta_a$','$\beta_{tp}$'}; % Parameter names
    else 
        error('Wrong system selected, choices: 1 = full, 2 = reduced');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  INPUT B) all other parameters, for Morris Method %%%%%%%%%%%%
elseif par_case == 2 % rest of the parameters
    if system == 1 % full system
        y = [gamma, Km, rM, rho0,alpha_m,...
           D_t,D_p,kappa_t,kappa_p,sph_t,sph_p,rt,rp, pos];
        Y_labels = {'$\gamma$', '$K_M$', '$r_M$', '$\rho_0$','$\alpha_m$',...
           '$D_t$','$D_p$','$\kappa_t$','$\kappa_p$','$s^{ph}_t$','$s^{ph}_p$',...
           '$r_t$','$r_p$','$x_1^{SF}$'};  % Parameter names
%         y = [gamma, Km, rM];
%         Y_labels = {'$\gamma$', '$K_M$', '$r_M$'};  % Parameter names
    elseif system == 2 % reduced system
        y = [beta_t,beta_p,D_t,D_p,kappa_t,kappa_p,sph_t,sph_p,rt,rp,...
            gamma, rM, Km,  pos];
        Y_labels = {'$D_t$','$D_p$','$\kappa_t$',...
            '$\kappa_p$','$s^{ph}_t$','$s^{ph}_p$',...
            '$r_t$','$r_p$','$\gamma$', '$r_M$', '$K_M$' ,'$x_1^{SF}$'};  % Parameter names
    else
        error('Wrong system selected, choices: 1 = full, 2 = reduced');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Parameter distributions and ranges of values
M = length(y);   % Number of uncertain parameters subject to SA
ymin = 0.01*y;
ymax = 100*y;
DistrFun = [repmat({'unif'}, M, 1)]; % Parameter distributions 
% To span many orders of magnitude this should be log-uniform, or log-normal
%%% Consider different distribution wider range for larger baseline values
% i.e. uniform distribution for rho0, kappa_t, kappa_p
% DistrFun(8) = {'unif'};  % rho0
% DistrFun(12) = {'unif'}; % kappa_t
% DistrFun(13) = {'unif'}; % kappa_p
% DistrFun(M) = {'unif'};  % pos must be unifromly distributed in [0,1]
DistrPar = cell(M,1);
for i=1:M
    DistrPar{i} = [ymin(i) ymax(i)];
end
if system == 2 && par_case == 1
    ymin(end) = 1.e-7;
    ymax(end) = 1.e-3;
    DistrPar{end} = [ymin(end) ymax(end)];  % cd must be unifromly distributed in [1e-9,1e-3]
end
if par_case == 2
    ymin(end) = 0;
    ymax(end) = 0.1;
    DistrPar{end} = [ymin(end) ymax(end)];  % pos must be uniformly distributed in [0,0.1]
end

%% %%%%%%%%%  SA 1) MORRIS METHOD, Elementary effects %%%%%%%%%%%%%%%%%%%%%
%                                                                            
%     - Output: active MMP concentration
%     - Input: (A) biochemical reaction rates in MMP-2 activation
%       and (B) all other parameters involved
%     - Distributions: log-normal (spanning several orders of magnitude)
%       or uniform (for larger parameters and those constrained in [0,1])
%     - Done for full system and reduced system to compare results
%     - One-at-a-time sampling, latin hypercube with radial sampling
%     - Elementary effects EE: mean calculated over abs(EE) to avoid
%       cancellation effects
%     - Convergence of the method checked via SAFE-inbuilt functions
%     - Uncertainty quantification via bootstrapping as provided in SAFE 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%  SAMPLING (One-at-a-time for Elementary effects)  %%%%%%%%%%%%


r = 500;   % Number of Elementary Effects
SampStrategy = 'rsu'; % Latin Hypercube (Alt: random uniform 'rsu') or lhs
% Ideal but not included in SAFE: Sobol sequence sampling
design_type = 'radial';

if compute_mod 
    %%% Sampling (One-at-a-time sampling for EE)
    Y = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);
    
    %%%%%%%%%%%%%%%%%%%%%%%%         EVALUATION        %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% First time
    Z = model_execution('PDE_solve',Y,test_case,system,par_case,L,Nx,x,dx,dim);
    %% Save
    save(strcat('../SA_saved/SA_EE_',num2str(system),'_r',num2str(r),'_',num2str(par_case),'.mat'),'r','Y','Z')
else
    load(['../SA_saved/SA_EE_',num2str(system),'_r',num2str(r),'_',num2str(par_case),'.mat'])
end

%%%%%%%%%%%%%%%   POST-PROCESSING: Elementary Effects Test  %%%%%%%%%%%%%%%%
%%


% load(strcat('SA_saved/SA_EE_',system,'_r',num2str(r),'.mat'),'r','Y','Z')

%%% EET: mean is calculated as the average of abs(EE)
[mi,sigma,EE] = EET_indices(r,ymin,ymax,Y,Z,design_type);
figure(1)
EET_plot(mi,sigma,Y_labels)
axis square

%%% Check convergence
rr = [ 1:1:r ] ;  
m_r = EET_convergence(EE,rr);   
figure
plot_convergence(m_r,rr*(M+1),[],[],[],'no of model evaluations','mean of EEs',Y_labels) 



%% Bootstrapping to derive confidence bounds

%%%%%%%%%%%%%%% Uncertainty quantification: Bootstrapping  %%%%%%%%%%%%%%%%

Nboot = 200; % Number of Bootstrap samples
[mi,sigma,EE,mi_sd,sigma_sd,mi_lb,sigma_lb,mi_ub,sigma_ub] = ...
    EET_indices(r,ymin,ymax,Y,Z,design_type,Nboot);



hfig = figure(77)
EET_plot(mi,sigma,Y_labels,mi_lb,mi_ub,sigma_lb,sigma_ub)
axis square


%% Repeat convergence analysis usign bootstrapping                                                                   
rr = [ 1:1:r ] ;  
[m_r,s_r,m_lb_r,m_ub_r] = EET_convergence(EE,rr,Nboot);                        
figure
plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],'no of model evaluations','mean of EEs',Y_labels)   



%% Mean and variance ( recall: mean is over abs(EE) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%       Boxplot      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = 'Helvetica' ; % font type of axes, labels, etc.
fs = 20 ; % font size of axes, labels, etc.
ms = 14 ; % marker size

hfig = figure(66)
data = [mi; sigma];
boxplot(data);
xticklabels(Y_labels)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
ylabel('Elementary Effects','FontSize',fs,'FontName',fn);
%title('Mean and Variance of Elementary Effects');

%% Visualise outliers (note: quartiles are computed with EE, not abs(EE) )
hfig = figure(55)
bp = boxplot(EE,'Symbol', 'o');
set(bp, 'MarkerSize', 7); % Adjust MarkerSize as needed
xticklabels(Y_labels)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
ylabel('Elementary Effects','FontSize',fs,'FontName',fn);

