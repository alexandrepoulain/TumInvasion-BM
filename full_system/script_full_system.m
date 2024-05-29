%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Script to simulate the full system. The conjunctive tissue is in
% 2D. the BM is the top boundary. 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
format long ;

addpath('../') % add path for function parameters

dim = 1; % the dimension of the conjunctive tissue (dimension fo the BM will be dim-1)
test = 3; % healthy = 1, tumor = 2, tumor+SF = 3

plot_SF = 1; % Parameter to maake aa plot of the effect of the SFs. 

%% Spatial discretization of conjunctive tissue
L = 0.1;   % Size of conjunctive tissue (length in dm)
Nx = 50; % Number of spatial grid points (in each direction if dim == 2)
dx = L/Nx; % the grid size
if dim == 1
    x = linspace(dx/2,L-dx/2,Nx);
    x_inter = linspace(0,L,Nx+1);
elseif dim == 2
    x = linspace(dx/2,L-dx/2,Nx);
    y = linspace(dx/2,L-dx/2,Nx);
    x_inter = linspace(0,L,Nx+1);
    y_inter = linspace(0,L,Nx+1);
    [X,Yx] = meshgrid(x,y);
else
    error("Wrong dimension")
end

%% Time parameters
Tf = 20;   % Final time (days)
Tp = 100; % Number of time points to store solution

%% Parameter values
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(test,L,dim,Nx,x);

% the source terms 
if dim == 1
    
    % Computation of the cell average of SF secretome 
    sf_secretome_t = @(x) (rt./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ));
    sf_secretome_p = @(x) (rp./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ));
    S_t = zeros(Nx, 1);
    S_p = zeros(Nx, 1);
    for ii=1:Nx
        S_t(ii) = sph_t + 1./dx*integral(sf_secretome_t,x_inter(ii),x_inter(ii+1)); 
        S_p(ii) = sph_p + 1./dx*integral(sf_secretome_p,x_inter(ii),x_inter(ii+1)); 
    end
elseif dim == 2  % pos is a matrix, nb line is the sumber of SFs, size 2 is 2  
    S_t = sph_t*ones(Nx*Nx, 1);
    S_p = sph_p*ones(Nx*Nx, 1);
    for ii = 1:size(posSF,1)
        %xSF =  reshape(X,[Nx^2,1]);
        xSF = X;
        %ySF =  reshape(Yx,[Nx^2,1]);
        ySF = Yx;
        
        sf_secretome_t = @(x,y) rt./(spread_secretum * sqrt(2*pi))*exp(...
                -(x-posSF(ii,1)).^2/(2*spread_secretum^2) -(y-posSF(ii,2)).^2/(2*spread_secretum^2) );
        sf_secretome_p = @(x,y) rp./(spread_secretum * sqrt(2*pi))*exp(...
                -(x-posSF(ii,1)).^2/(2*spread_secretum^2) -(y-posSF(ii,2)).^2/(2*spread_secretum^2) );
        
        

        for jj = 1:Nx % iterate on lines
            for rr = 1:Nx
                xmin = x_inter(rr);
                xmax = x_inter(rr+1);
                ymin = y_inter(jj);
                ymax = y_inter(jj+1);
                S_t((jj-1)*Nx+rr) = S_t((jj-1)*Nx+rr) + 1./(dx^2)*integral2(sf_secretome_t,xmin,xmax,ymin,ymax); 
                S_p((jj-1)*Nx+rr) = S_p((jj-1)*Nx+rr) +1./(dx^2)*integral2(sf_secretome_p,xmin,xmax,ymin,ymax); 
            end
        end
    end
end


%% optional: plot of SFs
if plot_SF && dim==2
    hfig = figure;  % save the figure handle in a variable
    S_t_plot = reshape(S_t,[Nx,Nx]);

    surf(X,Yx,S_t_plot,'LineWidth',3,'DisplayName','TIMP-2');
    view(2)
    shading interp
    clbr = colorbar
    xlabel('x axis')
    ylabel('y axis')

    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 0.8; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
    
    set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

end



%% Initial concentrations 
IC = initial_conditions(Nx,test,dim);

%% simulating the system
[t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
    = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
    width_BM );
