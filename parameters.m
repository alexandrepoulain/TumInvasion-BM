%%%% Baseline parameter set

function [gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,pos,rt,rp,spread_fib, alpha_rho_t, alpha_rho_p]...
    = parameters(test,L,dim,Nx,x)


% BM dynamics
gamma = 0.236; % s^-1
Km = 0.1357; % nM
rM = 6.18e-4; % nM/s
Mmax = 62.5e3; % nM

% Biochemical reaction rates
k0 = 2; % 1/(nM/s)
km0 = 1.e-2; % s^-1 
k1 = 2.71; % 1/(nM/s)
k2 = 0.14; % 1/(nM/s)
k3 = 0.02; % % s^(-1)
km1 = 1.e-4; % s^(-1)
km2 = 1.e-4; % s^(-1)

% Molecular decay rates
beta_t = 4.45e-5;
beta_ta = 0.035;
beta_m = 3.85e-2;
beta_d = 3.85e-2;
beta_p = 6.08e-4;
beta_a = 1.22e-3;
beta_tp = 3.60;

% tumor cells
if(test == 1)   
    rho0 = 0; 
elseif(test == 2| test == 3)
    if dim == 1
        rho0 = 1e+3; 
    elseif dim == 2
        rho0 = 1e+3*ones(Nx,1); % cells/dm^3
        rho0(x<L/4) = 0;
        rho0(x>3*L/4) = 0;
    end
else
    error("wrong test name")
end

% production of proMMP-2 and TIMP-2 by tumor cells
alpha_rho_t = 1e-5; 
alpha_rho_p = 1e-5; 

% Monomeric MMPs production
alpha_m = 5.*1e-4; % nmol/s/cell

% Diffusion rates in conjunctive tissue
D_t = 1.29e-6; % dm^2/s
D_p = 1.29e-6; 
% Diffusion rates in basement membrane
D_t_BM = D_t*1e-6; % dm^2/s
D_p_BM = D_p*1e-6;

% Transmission conditions
width_BM = 2e-6; % 200 nm
kappa_t = 3*D_t_BM/width_BM; 
kappa_p = 3*D_p_BM/width_BM;

% Physiological levels and sources
ct = 9.6;  
cp = 18.0;
sph_t = ct*beta_t;
sph_p = cp*beta_p;

% Senescent fibroblasts
if test == 3
    %rt = 2.22e-6; % TIMP-2 production by SF
    %rp = 2*4.14e-4; % proMMP-2 production by SF
    rt = 0.01*sph_t; % TIMP-2 production by SF
    rp = 0.08*sph_p; % proMMP-2 production by SF
else
    rt = 0;
    rp = 0;
end

nb_SF = 1; % number of SFs
spread_fib = 1e-2;

if dim == 1 % pos is a vector containing the position of each SF (x2)
    pos = zeros(nb_SF,1);
    for ii = 1:nb_SF
        pos(ii) = L-ii*0.005; % between 0 (furthest from BM) and 1 (at the BM)
    end
elseif dim == 2 % pos becomes a matrix, each line containing the position of one SF
    pos = zeros(nb_SF,2);
    for ii = 1:nb_SF
        pos(ii,1) = L-ii*0.039; % position in x1
        pos(ii,2) = L-ii*0.005; % position in x2
    end
end

end
