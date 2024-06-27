%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Main function to simulate the reduced system. Only in 1D. 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_BM,ct_conj, cp_conj, cp,ct,ca,M,cm,cd,c1,c2,c3,ctp,cta]...
    = reduced_system(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, width_BM, ...
    alpha_t, alpha_p)
%REDUCED_SYSTEM: Function to simulate the reduced system
    % Time parameters
    T = Tf*24*60*60; % Final time (in seconds)
    t_BM = linspace(0,T,Tp);
    
    % construct ODE system for diffusion equation in coonjunctive
    Lx = diffusion_matrix(Nx,dx, 2); % store laplacian matrix
    
    % Matrix for timp
    R_BC_timp = robin_bc(Nx,dx); % matrix to locate the nodes connected by transmission 
    R_BC_timp(1:Nx, 1:Nx+1) = kappa_t/dx * R_BC_timp(1:Nx, 1:Nx+1); % transmission boundary conditions for conjunctive
    R_BC_timp(Nx+1, 1: Nx+1) = kappa_t/width_BM*R_BC_timp(Nx+1, 1: Nx+1); % transmission for BM
    A_mat_timp = sparse(Nx+1,Nx+1); % initiate
    A_mat_timp(1:Nx,1:Nx) =  D_t*Lx; % the diffusion equation in conjunctive
    A_mat_timp = A_mat_timp + R_BC_timp; % add the transmission for both the BM and the conjunctive

    % Matrix for prommp
    R_BC_prommp = robin_bc(Nx,dx); % matrix to locate the nodes connected by transmission 
    R_BC_prommp(1:Nx, 1:Nx+1) = kappa_p/dx * R_BC_prommp(1:Nx, 1:Nx+1); % transmission boundary conditions for conjunctive
    R_BC_prommp(Nx+1, 1: Nx+1) = kappa_p/width_BM*R_BC_prommp(Nx+1, 1: Nx+1); % transmission for BM
    A_mat_prommp = sparse(Nx+1,Nx+1); % initiate
    A_mat_prommp(1:Nx,1:Nx) =  D_p*Lx; % the diffusion equation in conjunctive
    A_mat_prommp = A_mat_prommp + R_BC_prommp; % add the transmission for both the BM and the conjunctive
    
    Y0 = IC;
    
    %%% Solve ODE system until final time Tf
    options = odeset('Stats','on', 'RelTol',1e-8,'AbsTol',1e-10);
    [t,Y] = ode15s(@(t,Y) ODEhealthyredsysfun(t,Y, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, A_mat_timp, A_mat_prommp, Nx, alpha_t,alpha_p ), ...
        t_BM, Y0, options);
        
        
    ct_conj = Y(:,1:Nx);
    cp_conj = Y(:,Nx+1:2*Nx);
    
    M = Y(:,2*Nx+1);
    cm = Y(:,2*Nx+2);
    cd = Y(:,2*Nx+3);
    cp = Y(:,2*Nx+4);
    ct = Y(:,2*Nx+5);
    ca = Y(:,2*Nx+6);
    
    c1 = Y(:,2*Nx+7);
    c2 = Y(:,2*Nx+8);
    c3 = Y(:,2*Nx+9);

    ctp = Y(:,2*Nx+10);
    cta = Y(:,2*Nx+11);

    t_BM = t;

end 

function [dYdt] = ODEhealthyredsysfun(t, Y, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, A_mat_timp, A_mat_prommp, Nx ...
    , alpha_t,alpha_p)

    dYdt = zeros(Nx*2+11,1);

    ct_conj = Y(1:Nx);
    cp_conj = Y(Nx+1:2*Nx);
    
    M = Y(2*Nx+1);
    cm = Y(2*Nx+2);
    cd = Y(2*Nx+3);
    cp = Y(2*Nx+4);
    ct = Y(2*Nx+5);
    ca = Y(2*Nx+6);
    
    c1 = Y(2*Nx+7);
    c2 = Y(2*Nx+8);
    c3 = Y(2*Nx+9);

    ctp = Y(2*Nx+10);
    cta = Y(2*Nx+11);
    
    
    ct_tot = [ct_conj;ct];
    cp_tot = [cp_conj;cp];
    
    dct_conj = zeros(Nx,1); 
    dcp_conj = zeros(Nx,1);
    % the diffusion and transmission
    dct_buff= A_mat_timp*ct_tot;
    dcp_buff = A_mat_prommp*cp_tot;
    
    dct_conj = dct_buff(1:Nx);
    dct = dct_buff(Nx+1);
    
    dcp_conj = dcp_buff(1:Nx);
    dcp = dcp_buff(Nx+1);

    % the reactions
    dct_conj = dct_conj + S_t - beta_t*ct_conj;
    dcp_conj = dcp_conj + S_p - beta_p*cp_conj;
% 
%     dct_conj = zeros(Nx,1);
%     dcp_conj = zeros(Nx,1);
    
    
    c1 = k1*cd/km1*ct; 
    c2 = k1*k2*cd*ct*cp/(km1*(km2+k3));
    c3 = (k1/km1)^2*cd*ct^2;
    ctp = k2*ct*cp/beta_tp;
    cta = k2*ca*ct/beta_ta;
    
    dM = - gamma * ca * M /(Km + M) + rM*(1-M/Mmax);
    dcm = 0;
    dcd = 0;

    dcp = dcp  - k2 *cp.*ct*(1 + k1*k3/(km1*(km2+k3))*cd) ; % + rho0*alpha_p
    dct = dct  - ct.*k2*(cp+ ca) ; %  + rho0*alpha_t
    dca = (k3*c2 - k2*ct.*ca) ;
    dc1 = 0;
    dc2 = 0;
    dc3 = 0;
    dctp = 0;
    dcta = 0;
 
    dYdt = [dct_conj; dcp_conj; dM;dcm; dcd; dcp;dct;dca; dc1;dc2;dc3;dctp;dcta];

end
