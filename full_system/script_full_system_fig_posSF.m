%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Script to simulate the full system and change the position of the SF.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
format long ;

addpath('../') % add path for function parameters

dim = 1; % the dimension of the conjunctive tissue (dimension fo the BM will be dim-1)
test = 3; % % healthy = 1, tumor = 2, tumor+SF = 3
plot_SF = 0;

if dim == 1
    addpath('./1D/') 
elseif dim == 2
     addpath('./2D/') 
else
    error('Wrong dimension, select 1 or 2')
end

% Spatial discretization of conjunctive tissue
L = 0.1;   % Size of conjunctive tissue (length in dm)
Nx = 40; % Number of spatial grid points (in each direction if dim == 2)
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

% Time parameters
Tf = 2;   % Final time (days)
Tp = 100; % Number of time points to store solution

% Parameter values
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(test,L,dim,Nx,x);

% Variable SF position along the line and sph_t
sph_t_vec = logspace(-10,-3,100);
sph_p_vec = logspace(-7,0,100);

%rp_vec = logspace(-9,-2,100);

%posSF_vec = linspace(0,L,100);
rt = 0
posSF = 0.005;
length(sph_t_vec)

for iii=1:length(sph_t_vec)
    iii
    for jjj=1:length(sph_p_vec)
        sph_t = sph_t_vec(iii);
        %posSF = posSF_vec(jjj);
        sph_p = sph_p_vec(jjj);
        %rp = rp_vec(jjj);
        % the source terms 
        if dim == 1
            S_t = (sph_t + rt./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ))'; 
            S_p = (sph_p + rp./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ))'; 
        elseif dim == 2  % pos is a matrix, nb line is the sumber of SFs, size 2 is 2  
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
        
        
        % optional: plot of SFs
        if plot_SF && dim==2
            hfig = figure;  % save the figure handle in a variable
            S_t_plot = reshape(S_t,[Nx,Nx]);
        
            surf(X,Yx,S_t_plot,'LineWidth',3,'DisplayName','TIMP-2');
            view(2)
            shading interp
            clbr = colorbar
            xlabel('x axis')
            ylabel('y axis')
        
            %fname = strcat(test , '_concentrations_BM.eps');
            %fname = strcat(plotpath,fname)
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig,'-property','Box'),'Box','off') % optional
            set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig,'Position');
            set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            %legend('Location', 'east');
        end
        
        
        
        % Initial concentrations 
        IC = initial_conditions(Nx,test,dim);
        
        % Test for stability
        test_stability(kappa_p,kappa_t,k1,k2,k3,km1,km2,Km,gamma,rM,Mmax,...
            IC(2*Nx+5),IC(2*Nx+4),IC(2*Nx+6),IC(2*Nx+3), S_t(end), S_p(end),...
            rho0, alpha_t,alpha_p,width_BM, beta_t, beta_p)
        
        % simulating the system
        [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
            = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
            rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
            beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
            width_BM );

        % Saving the last BM density in matrix
        BM_dens_mat(iii,jjj) = M_BM(end);
        dict_BM(((iii-1)*length(sph_t_vec)+1)+jjj,:) = [sph_t,sph_p, M_BM(end)];
    end
end
%%
[X,Y] = meshgrid(sph_p_vec,sph_t_vec)

if test == 1
    test_case_str = "healthy";
elseif test ==2
    test_case_str = "tumor";
elseif test == 3
    test_case_str = "tumor+SF";
end
plotpath = '../results/full/1D/'
fname = strcat(test_case_str , '_MBdens_fnof_spht_sphp.eps');
fname = strcat(plotpath,fname)
hfig = figure;  % save the figure handle in a variable

surf(X,Y,BM_dens_mat,'DisplayName','BM density'); 
view(2)
shading interp
clbr = colorbar
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

title('$M$ (nM)')
%ylabel(clbr,,'Rotation',270)
grid on;
ylabel('$s_{ph}^t$ (nM/s)')
xlabel('$s_{ph}^p$ (nM/s)')
xlim([1e-7, 1])
ylim([1e-10, 1e-3])

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
clbr.Label.Interpreter = 'latex'

set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
%legend('Location', 'east');
saveas(hfig,fname, 'epsc')
