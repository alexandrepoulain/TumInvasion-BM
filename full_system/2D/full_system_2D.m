function [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
    = full_system_2D(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, width_BM )
%FULL_SYSTEM_2D  Function to simulate the full system in 2D 
    
    % Time parameters
    T = Tf*24*60*60; % Final time (in seconds)
    t_BM = linspace(0,T,Tp);
    
    % construct system in matrix form
    Lx = diffusion_matrix(Nx,dx); % store laplacian matrix
    R_BC_t = robin_bc(Nx,dx); % matrix to locate the nodes connected by transmission 
    R_BC_t(1:Nx^2, 1:Nx^2+Nx) = kappa_t/dx * R_BC_t(1:Nx^2, 1:Nx^2+Nx); % transmission boundary conditions for conjunctive
    R_BC_t(Nx^2+1:end, 1: end) = kappa_t/width_BM*R_BC_t(Nx^2+1:end, 1: end); % transmission for BM
    R_BC_p = robin_bc(Nx,dx); % matrix to locate the nodes connected by transmission 
    R_BC_p(1:Nx^2, 1:Nx^2+Nx) = kappa_p/dx * R_BC_p(1:Nx^2, 1:Nx^2+Nx); % transmission boundary conditions for conjunctive
    R_BC_p(Nx^2+1:end, 1: end) = kappa_p/width_BM*R_BC_p(Nx^2+1:end, 1: end); % transmission for BM

    % create matrix associated to the system
    A_t = sparse(Nx^2+Nx,Nx^2+Nx); % initiate
    A_t(1:Nx^2,1:Nx^2) = D_t*Lx; % the diffusion equation in conjunctive
    A_t = A_t + R_BC_t; % add the transmission for both the BM and the conjunctive
    A_p = sparse(Nx^2+Nx,Nx^2+Nx); % initiate
    A_p(1:Nx^2,1:Nx^2) = D_p*Lx; % the diffusion equation in conjunctive
    A_p = A_p + R_BC_p; % add the transmission for both the BM and the conjunctive

    Y0 = IC;
    options = odeset('Stats','on');
    [t,Y] = ode15s(@(t,Y) ODEfullsysfun_2D(t,Y, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, A_t, A_p, Nx ), ...
        t_BM, Y0, options);

    for ii = 1:length(t)

        ct_conj(ii,:) = Y(ii,1:Nx^2);
        cp_conj(ii,:) = Y(ii,Nx^2+1:2*Nx^2);
        
        M(ii,:) = Y(ii,2*Nx^2+1: 2*Nx^2+Nx);
        cm_BM(ii,:) = Y(ii,2*Nx^2+Nx+1: 2*Nx^2+2*Nx);
        cd_BM(ii,:) = Y(ii,2*Nx^2+2*Nx+1: 2*Nx^2+3*Nx);
        cp_BM(ii,:) = Y(ii,2*Nx^2+3*Nx+1: 2*Nx^2+4*Nx);
        ct_BM(ii,:) = Y(ii,2*Nx^2+4*Nx+1: 2*Nx^2+5*Nx);
        ca_BM(ii,:) = Y(ii,2*Nx^2+5*Nx+1: 2*Nx^2+6*Nx);
        
        c1_BM(ii,:) = Y(ii,2*Nx^2+6*Nx+1: 2*Nx^2+7*Nx);
        c2_BM(ii,:) = Y(ii,2*Nx^2+7*Nx+1: 2*Nx^2+8*Nx);
        c3_BM(ii,:) = Y(ii,2*Nx^2+8*Nx+1: 2*Nx^2+9*Nx);
    
        ctp_BM(ii,:) = Y(ii,2*Nx^2+9*Nx +1 : 2*Nx^2+10*Nx);
        cta_BM(ii,:) = Y(ii,2*Nx^2+10*Nx+1: 2*Nx^2+11*Nx);
    end
    
    t_BM = t;

end


function [dYdt] = ODEfullsysfun_2D(t, Y, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, A_mat_timp, A_mat_prommp, Nx )
    % This function describes the ODE system associated to diffusion and 
    % transmission to BM
    t
    dYdt = zeros(Nx^2*2+Nx*11,1);

    ct_conj = Y(1:Nx^2);
    cp_conj = Y(Nx^2+1:2*Nx^2);
    
    M = Y(2*Nx^2+1: 2*Nx^2+Nx);
    cm = Y(2*Nx^2+Nx+1: 2*Nx^2+2*Nx);
    cd = Y(2*Nx^2+2*Nx+1: 2*Nx^2+3*Nx);
    cp = Y(2*Nx^2+3*Nx+1: 2*Nx^2+4*Nx);
    ct = Y(2*Nx^2+4*Nx+1: 2*Nx^2+5*Nx);
    ca = Y(2*Nx^2+5*Nx+1: 2*Nx^2+6*Nx);
    
    c1 = Y(2*Nx^2+6*Nx+1: 2*Nx^2+7*Nx);
    c2 = Y(2*Nx^2+7*Nx+1: 2*Nx^2+8*Nx);
    c3 = Y(2*Nx^2+8*Nx+1: 2*Nx^2+9*Nx);

    ctp = Y(2*Nx^2+9*Nx +1 : 2*Nx^2+10*Nx);
    cta = Y(2*Nx^2+10*Nx+1: 2*Nx^2+11*Nx);
   
    ct_tot = [ct_conj;ct];
    cp_tot = [cp_conj;cp];
    
    % init
    dct_conj = zeros(Nx^2,1); 
    dcp_conj = zeros(Nx^2,1);
    % the diffusion and transmission
    dct_buff= A_mat_timp*ct_tot;
    dcp_buff = A_mat_prommp*cp_tot;
    
    dct_conj = dct_buff(1:Nx^2);
    dct = dct_buff(Nx^2+1:Nx^2+Nx);
    
    dcp_conj = dcp_buff(1:Nx^2);
    dcp = dcp_buff(Nx^2+1:Nx^2+Nx);

    % the reactions
    dct_conj = dct_conj + S_t - beta_t*ct_conj;
    dcp_conj = dcp_conj + S_p - beta_p*cp_conj;
    
    dM = - gamma * ca .* M ./(Km + M) + rM.*(1-M/Mmax);
    dcm = alpha_m*rho0 - k0.*power(cm,2) + km0*cd - beta_m * cm;
    dcd = k0.*power(cm,2) - km0*cd - k1*cd.*ct + km1*c1 - beta_d*cd;
    dcp = dcp - k2 * cp.*(ct+c1) + km2*(ctp + c2) - beta_p*cp;
    dct = dct - ct.*(k1*(cd+c1) +k2*(cp+ca)) + km1*(c1+c3) + km2*(ctp + cta) -beta_t*ct;
    dca = k3*c2 - k2*ct.*ca + km2*cta - beta_a*ca;
    dc1 = k1*cd.*ct - c1.*(k1*ct +k2*cp) + km1*(c3-c1) + (km2+k3)*c2- beta_d*c1;
    dc2 = k2*c1.*cp - (km2+k3)*c2- beta_d*c2;
    dc3 = k1*c1.*ct - km1*c3 - beta_d*c3;
    dctp = k2*ct.*cp - km2*ctp - beta_tp *ctp;
    dcta = k2*ct.*ca - km2*cta - beta_ta*cta;
 
    dYdt = [dct_conj; dcp_conj; dM;dcm; dcd; dcp;dct;dca; dc1;dc2;dc3;dctp;dcta];
end
