%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% function to simulate the full system. 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t_BM,ct_conj, cp_conj, cp,ct,ca,M,cm,cd,c1,c2,c3,ctp,cta] =  ... 
    full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
    width_BM )
%FULL_SYSTEM This function calls the correct function depending on the
%dimension of the problem
if dim == 1
    [t_BM,ct_conj, cp_conj, cp,ct,ca,M,cm,cd,c1,c2,c3,ctp,cta]...
    = full_system_1D(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_t, D_p, ...
    width_BM );
elseif dim == 2
    [t_BM,ct_conj, cp_conj, cp,ct,ca,M,cm,cd,c1,c2,c3,ctp,cta]...
    = full_system_2D(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
    rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
    beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_t, D_p, ...
    width_BM );
else
    error("Cannot continue... dimension incorrect")
end

end

