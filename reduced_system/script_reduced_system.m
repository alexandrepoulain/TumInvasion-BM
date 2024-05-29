%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to simulate the full system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Thu Oct 23 2023
% 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

format long ;

addpath('../')

dim = 1; % the dimension of the conjunctive tissue (if healthy case, always 1d)
test = 2; % test case "healthy" = 1 , "tumor" = 2, 'tumor+ SF' = 3

% Spatial discretization of conjunctive tissue
L = 0.1;   % Size of conjunctive tissue (length in dm)
Nx = 100; % Number of spatial grid points (1D cross section conjuntive tissue)
dx = L/Nx; 
x = linspace(0,L,Nx);
% Time parameters
Tf = 1;   % Final time (days)
Tp = 100; % Number of time points to store solution

% Parameter values
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,pos,rt,rp,spread_secretum,alpha_t, alpha_p] ...
    = parameters(test,L,dim,Nx,x);

% the source terms 
S_t = (sph_t + rt./(spread_secretum * sqrt(2*pi))*exp(-(x-pos).^2/(2*spread_secretum^2) ))'; 
S_p = (sph_p + rp./(spread_secretum * sqrt(2*pi))*exp(-(x-pos).^2/(2*spread_secretum^2) ))'; 

% Initial concentrations 
IC = initial_conditions(Nx,test);

% Test for stability
test_stability(kappa_p,kappa_t,k1,k2,k3,km1,km2,Km,gamma,rM,Mmax,...
    IC(2*Nx+5),IC(2*Nx+4),IC(2*Nx+6),IC(2*Nx+3), S_t(end), S_p(end),...
    rho0, alpha_t,alpha_p,width_BM)

% simulating the system
[t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
    = reduced_system(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_t, D_p, width_BM,alpha_t, alpha_p );

c1_BM = k1*cd_BM/km1.*ct_BM; 
c2_BM = k1*k2*cd_BM.*ct_BM.*cp_BM/(km1*(km2+k3));
c3_BM = (k1/km1)^2*cd_BM.*ct_BM.^2;


